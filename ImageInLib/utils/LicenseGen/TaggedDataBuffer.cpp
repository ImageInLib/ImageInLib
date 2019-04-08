// TaggedDataBuffer.cpp: implementation of the CTaggedDataBuffer class.
//
//////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "TaggedDataBuffer.h"

//====================================================================
//
// Element header:
//
//  bit 8  7  6  5  4  3  2  1
//     (X)(T  T)(L  L  L  L  L) <tag> [<length>] [<data>]
//
//      X - element type: 0 = value
//                        1 = sequence
//
//     TT - tag length:   00 - 1 byte
//                        01 - 2 bytes
//                        10 - 3 bytes
//                        11 - 4 bytes
//
//  LLLLL - data lenght   0-27 - implict length 1-28 bytes
//                        28   - 1 byte
//                        29   - 2 bytes
//                        30   - 3 bytes
//                        31   - 4 bytes
//
//
// Sequence item encoding:
//
//  bit 8  7  6  5  4  3  2  1
//     (L  L  L  L  L  L  L  L) [<length>] [<data>]
//
//  LLLLLLLL - 0-252 - implicit length 0-252 bytes
//             253   - 2 bytes
//             254   - 3 bytes
//             255   - 4 bytes
//
//====================================================================


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTaggedDataBuffer buf;


CTaggedDataBuffer::CTaggedDataBuffer()
{
	Init();
}

CTaggedDataBuffer::~CTaggedDataBuffer()
{
}


void CTaggedDataBuffer::Clear()
{
	Init();
}

void CTaggedDataBuffer::Init()
{
	m_RootElement.items.clear();
	m_RootElement.dataAllocated = false;
	m_RootElement.dwLength = 0;
	m_RootElement.dwTag = 0;
	m_RootElement.dwHeaderSize = 0;
	m_RootElement.pData = NULL;
	m_RootElement.pParent = NULL;
	m_RootElement.fragmented = false;
	m_RootElement.items.push_back(new SequenceItem);
	m_RootElement.items[0]->dwHeaderSize = 0;
	m_RootElement.items[0]->dwLength = 0;
	m_RootElement.items[0]->pData = NULL;
	m_RootElement.items[0]->itemIdx = 0;
	m_RootElement.items[0]->pParent = &m_RootElement;
	m_RootElement.items[0]->fragmented = false;
}

bool CTaggedDataBuffer::Parse(const BYTE *pBuffer, DWORD dwBufferSize)
{
	bool bRet;

	if (m_RootElement.pData != NULL) {
		delete m_RootElement.pData;
		m_RootElement.pData = NULL;
	}
	m_RootElement.pData = new BYTE[dwBufferSize];
	if (m_RootElement.pData == NULL)
		return false;
	memmove(m_RootElement.pData, pBuffer, dwBufferSize);
	m_RootElement.dwLength = dwBufferSize;
	m_RootElement.dataAllocated = true;

	bRet = ParseSequence(&m_RootElement);

	if (!bRet)
		Init();

	return bRet;
}


bool CTaggedDataBuffer::ParseSequence(Element *pSequence)
{
	stack<ParserState> stateStack;
	ParserState parser;
	ElementType elementType;
	DWORD headerSize; 
	DWORD elementTag; 
	DWORD elementLength;
	BYTE *pElementData;

	pSequence->fragmented = false;
	pSequence->items.clear();
	if (pSequence->dwLength == 0)
		return true;

	if (pSequence != &m_RootElement) {
		parser.seq = pSequence;
		parser.buf.dwRemaining = 0;
		parser.seqBuf.dwRemaining = pSequence->dwLength;
		parser.seqBuf.pCursor = pSequence->pData;
		parser.lastTag = 0;
	} else {
		pSequence->items.push_back(new SequenceItem);
		pSequence->items[0]->dwHeaderSize = 0;
		pSequence->items[0]->dwLength = 0;
		pSequence->items[0]->pData = pSequence->pData;
		pSequence->items[0]->itemIdx = 0;
		pSequence->items[0]->pParent = pSequence;
		pSequence->items[0]->fragmented = false;
		parser.buf.pCursor = pSequence->pData;
		parser.buf.dwRemaining = pSequence->dwLength;
		parser.lastTag = 0;
		parser.seq = pSequence;
		parser.elements = &pSequence->items[0]->elements;
	}

	for (;;) {
		if (GetAvailable(parser.buf) == 0) {
			if (parser.seq != &m_RootElement && GetAvailable(parser.seqBuf) != 0) {
				if (!ParseSequenceItem(parser.seqBuf, headerSize, elementLength, pElementData))
					return false;
				SequenceItem *item = new SequenceItem;
				item->dwHeaderSize = headerSize;
				item->dwLength = elementLength;
				item->pData = pElementData;
				item->itemIdx = parser.seq->items.size();
				item->fragmented = false;
				item->pParent = parser.seq;
				parser.seq->items.push_back(item);
				parser.elements = &item->elements;
				parser.buf.dwRemaining = elementLength;
				parser.buf.pCursor = pElementData;
				parser.lastTag = 0;
				continue;
			}
			if (stateStack.size() == 0)
				break;
			parser = stateStack.top();
			stateStack.pop();
			continue;
		}

		if (!ParseElement(parser.buf, elementType, headerSize, elementTag, elementLength, pElementData))
			return false;

		parser.lastTag += elementTag;
		elementTag = parser.lastTag;
		if (parser.elements->find(elementTag) != parser.elements->end())
			return false;
		Element &element = (*parser.elements)[elementTag];

		element.type = elementType;
		element.dwTag = elementTag;
		element.dwHeaderSize = headerSize;
		element.dwLength = elementLength;
		element.pData = pElementData;
		element.dataAllocated = false;
		element.fragmented = false;
		element.pParent = parser.seq;
		element.itemIdx = parser.seq->items.size() - 1;

		if (elementType == etSequence) {
			stateStack.push(parser);
			parser.seq = &element;
			parser.buf.dwRemaining = 0;
			parser.seqBuf.dwRemaining = elementLength;
			parser.seqBuf.pCursor = pElementData;
			parser.lastTag = 0;
		}
	}

	return true;
}


bool CTaggedDataBuffer::Write(BYTE *pBuffer, DWORD *pdwBufferSize)
{
	if (!DefragmentBuffer())
		return false;
	*pdwBufferSize = m_RootElement.items[0]->dwLength;
	if (pBuffer && m_RootElement.pData)
		memmove(pBuffer, m_RootElement.items[0]->pData, m_RootElement.items[0]->dwLength);
	return true;
}


bool CTaggedDataBuffer::ParseElement(
	BufferState &bufState, 
	ElementType &elementType,
	DWORD &headerSize, 
	DWORD &elementTag, 
	DWORD &elementLength, 
	BYTE* &pElementData)
{
	BYTE hdr;
	BYTE *pData;
	DWORD length;

	if (!GetByte(bufState, hdr))
		return false;

	headerSize = 1;

	if (hdr & 0x80) 
		elementType = etSequence;
	else
		elementType = etValue;

	length = ((hdr & 0x60) >> 5) + 1;
	headerSize += length;
	if (!GetBytes(bufState, length, pData))
		return false;
	if (length == 1)
		elementTag = *pData;
	else if (length == 2)
		elementTag = pData[0] + (pData[1] << 8);
	else if (length == 3)
		elementTag = pData[0] + (pData[1] << 8) + (pData[2] << 16);
	else
		elementTag = pData[0] + (pData[1] << 8) + (pData[2] << 16) + (pData[3] << 24);

	elementLength = hdr & 0x1f;
	if (elementLength > 27) {
		if (!GetBytes(bufState, elementLength - 27, pData))
			return false;
		headerSize += elementLength - 27;
		if (elementLength == 28)
			elementLength = *pData;
		else if (elementLength == 29)
			elementLength = pData[0] + (pData[1] << 8);
		else if (elementLength == 30)
			elementLength = pData[0] + (pData[1] << 8) + (pData[2] << 16);
		else
			elementLength = pData[0] + (pData[1] << 8) + (pData[2] << 16) + (pData[3] << 24);
	} else {
		elementLength++;
	}

	if (!GetBytes(bufState, elementLength, pElementData))
		return false;

	return true;
}


bool CTaggedDataBuffer::ParseSequenceItem(
	BufferState &bufState, 
	DWORD &headerSize, 
	DWORD &itemLength, 
	BYTE* &pItemData)
{
	BYTE hdr;
	BYTE *pData;


	if (!GetByte(bufState, hdr))
		return false;

	headerSize = 1;

	itemLength = hdr;
	if (itemLength > 252) {
		if (!GetBytes(bufState, itemLength - 251, pData))
			return false;
		headerSize += itemLength - 251;
		if (itemLength == 253)
			itemLength = pData[0] + (pData[1] << 8);
		else if (itemLength == 254)
			itemLength = pData[0] + (pData[1] << 8) + (pData[2] << 16);
		else
			itemLength = pData[0] + (pData[1] << 8) + (pData[2] << 16) + (pData[3] << 24);
	}

	if (!GetBytes(bufState, itemLength, pItemData))
		return false;

	return true;	
}


HTDBITEM CTaggedDataBuffer::GetItem(HTDBELEMENT hElement, DWORD dwIndex)
{
	Element *pElement = dynamic_cast<Element*>(hElement);

	if (pElement == NULL)
		return NULL;

	if (pElement->type != etSequence)
		return NULL;

	if (pElement->items.size() <= dwIndex)
		return NULL;

	return pElement->items[dwIndex];
}


HTDBITEM CTaggedDataBuffer::GetNextItem(HTDBITEM hItem)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);

	if (pItem == NULL)
		return NULL;

	if (pItem->itemIdx == pItem->pParent->items.size()-1)
		return NULL;

	return pItem->pParent->items[pItem->itemIdx+1];
}


DWORD CTaggedDataBuffer::GetItemsCount(HTDBELEMENT hElement)
{
	Element *pElement = dynamic_cast<Element*>(hElement);

	if (pElement == NULL)
		return (DWORD)-1;

	if (pElement->type != etSequence)
		return (DWORD)-1;

	return pElement->items.size();
}


HTDBELEMENT CTaggedDataBuffer::GetFirstElement(HTDBITEM hItem)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);

	if (hItem == NULL) {
		if (m_RootElement.items[0]->elements.size() == 0)
			return NULL;
		return &m_RootElement.items[0]->elements.begin()->second;
	}

	if (pItem == NULL)
		return NULL;

	if (pItem->elements.size() == 0)
		return NULL;
	return &pItem->elements.begin()->second;
}


HTDBELEMENT CTaggedDataBuffer::GetNextElement(HTDBELEMENT hElement)
{
	Element *pElement = dynamic_cast<Element*>(hElement);

	if (pElement == NULL)
		return NULL;

	ElementMap &elements = pElement->pParent->items[pElement->itemIdx]->elements;
	ElementMap::iterator it = elements.find(pElement->dwTag);
	if (it == elements.end())
		return NULL;

	return &(++it)->second;
}


HTDBELEMENT CTaggedDataBuffer::FindElement(HTDBITEM hItem, DWORD dwTag)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);
	ElementMap *pElements;

	if (hItem == NULL)
		pElements = &m_RootElement.items[0]->elements;
	else
		pElements = &pItem->elements;

	ElementMap::iterator it = pElements->find(dwTag);
	if (it == pElements->end())
		return NULL;
	return &it->second;
}


bool CTaggedDataBuffer::GetItemInfo(HTDBITEM hItem, DWORD* pdwTag, DWORD* pdwIndex, DWORD* pdwLength, BYTE** ppData)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);

	if (hItem == NULL)
		pItem = m_RootElement.items[0];

	if (pItem == NULL)
		return false;

	if (pItem->fragmented) {
		if (!DefragmentBuffer())
			return false;
	}

	if (pdwTag)
		*pdwTag = pItem->pParent->dwTag;
	if (pdwIndex)
		*pdwIndex = pItem->itemIdx;
	if (pdwLength)
		*pdwLength = pItem->dwLength;
	if (ppData)
		*ppData = pItem->pData;

	return true;
}


bool CTaggedDataBuffer::GetElementInfo(HTDBELEMENT hElement, ElementType* pType, DWORD* pdwTag, DWORD* pdwLength, BYTE** ppData)
{
	Element *pElement = dynamic_cast<Element*>(hElement);

	if (pElement == NULL)
		return false;

	if (pElement->fragmented) {
		if (!DefragmentBuffer())
			return false;
	}

	if (pType)
		*pType = pElement->type;
	if (pdwTag)
		*pdwTag = pElement->dwTag;
	if (pdwLength)
		*pdwLength = pElement->dwLength;
	if (ppData)
		*ppData = pElement->pData;

	return true;
}


bool CTaggedDataBuffer::SetElementValue(HTDBITEM hItem, ElementType nType, DWORD dwTag, DWORD dwLength, const BYTE* pData)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);
	Element *pElement;
	BYTE *pDataBuffer = NULL;

	if (hItem == NULL)
		pItem = m_RootElement.items[0];

	if (pItem == NULL)
		return false;

	pDataBuffer = new BYTE[dwLength];
	if (pDataBuffer == NULL)
		return false;
	memmove(pDataBuffer, pData, dwLength);

	ElementMap::iterator it = pItem->elements.find(dwTag);
	if (it != pItem->elements.end()) {
		pElement = &it->second;
		pElement->items.clear();
		pElement->fragmented = false;
	} else {
		pElement = &pItem->elements[dwTag];
		pElement->dataAllocated = false;
		pElement->fragmented = false;
		pElement->dwLength = 0;
		pElement->dwTag = dwTag;
		pElement->dwHeaderSize = 0;
		pElement->itemIdx = pItem->itemIdx;
		pElement->pData = NULL;
		pElement->pParent = pItem->pParent;
	}

	if (pElement->dataAllocated)
		delete pElement->pData;
	pElement->dataAllocated = true;
	pElement->pData = pDataBuffer;
	pElement->dwLength = dwLength;
	pElement->type = nType;

	if (nType == etSequence && !ParseSequence(pElement))
		return false;

	for (;;) {
		if (pItem->pParent->fragmented)
			break;
		pItem->fragmented = true;
		pItem->pParent->fragmented = true;
		if (pItem->pParent == &m_RootElement)
			break;
		pItem = pItem->pParent->pParent->items[pItem->pParent->itemIdx];
	}

	return true;
}


HTDBELEMENT CTaggedDataBuffer::CreateSequence(HTDBITEM hItem, DWORD dwTag)
{
	SequenceItem *pItem = dynamic_cast<SequenceItem*>(hItem);
	Element *pElement;
	BYTE *pDataBuffer = NULL;

	if (hItem == NULL)
		pItem = m_RootElement.items[0];

	if (pItem == NULL)
		return false;

	ElementMap::iterator it = pItem->elements.find(dwTag);
	if (it != pItem->elements.end()) {
		return &it->second;
	} else {
		pElement = &pItem->elements[dwTag];
		pElement->dataAllocated = false;
		pElement->fragmented = true;
		pElement->dwLength = 0;
		pElement->dwTag = dwTag;
		pElement->dwHeaderSize = 0;
		pElement->itemIdx = pItem->itemIdx;
		pElement->pData = NULL;
		pElement->pParent = pItem->pParent;
		pElement->type = etSequence;
	}

	for (;;) {
		if (pItem->pParent->fragmented)
			break;
		pItem->fragmented = true;
		pItem->pParent->fragmented = true;
		if (pItem->pParent == &m_RootElement)
			break;
		pItem = pItem->pParent->pParent->items[pItem->pParent->itemIdx];
	}

	return pElement;
}


HTDBITEM CTaggedDataBuffer::CreateSequenceItem(HTDBELEMENT hSequence)
{
	Element *pSequence = dynamic_cast<Element*>(hSequence);

	if (pSequence == NULL)
		return NULL;
	if (pSequence->type != etSequence)
		return NULL;

	SequenceItem *item = new SequenceItem;
	item->fragmented = true;
	item->dwLength = 0;
	item->pData = NULL;
	item->itemIdx = pSequence->items.size();
	item->pParent = pSequence;
	pSequence->items.push_back(item);
	
	SequenceItem *pItem = pSequence->items.back();
	for (;;) {
		if (pItem->pParent->fragmented)
			break;
		pItem->fragmented = true;
		pItem->pParent->fragmented = true;
		if (pItem->pParent == &m_RootElement)
			break;
		pItem = pItem->pParent->pParent->items[pItem->pParent->itemIdx];
	}

	return pItem;
}


bool CTaggedDataBuffer::DefragmentBuffer()
{
	stack<DefragState> stateStack;
	DefragState state;
	BufferState buf;
	DWORD dwDataBufferLength;
	BYTE *pDataBuffer;

	if (!m_RootElement.fragmented)
		return true;

	if (m_RootElement.items.size() == 0) {
		if (m_RootElement.dataAllocated)
			delete m_RootElement.pData;
		m_RootElement.pData = NULL;
		m_RootElement.dwLength = 0;
		return true;
	}

	state.seq = &m_RootElement;
	state.itemIt = m_RootElement.items.begin();
	state.elemIt = (*state.itemIt)->elements.begin();
	state.lastTag = 0;

	// update lengths
	state.seq->dwLength = 0;
	if (state.itemIt != state.seq->items.end())
		(*state.itemIt)->dwLength = 0;
	for (;;) {
		SequenceItem *pItem = *state.itemIt;
		if (state.elemIt != pItem->elements.end()) {
			Element &element = state.elemIt->second;
			if (element.type == etSequence && element.fragmented) {
				if (element.items.size() == 0) {
					element.dwLength = 0;
					element.dwHeaderSize = CalculateElementHeaderSize(element.dwTag, state.lastTag, element.dwLength);
					pItem->dwLength += element.dwHeaderSize + element.dwLength;
					state.lastTag = state.elemIt->second.dwTag;
					state.elemIt++;
				} else {
					stateStack.push(state);
					state.seq = &element;
					state.seq->dwLength = 0;
					state.itemIt = state.seq->items.begin();
					state.elemIt = (*state.itemIt)->elements.begin();
					state.lastTag = 0;
				}
			} else {
				element.dwHeaderSize = CalculateElementHeaderSize(element.dwTag, state.lastTag, element.dwLength);
				pItem->dwLength += element.dwHeaderSize + element.dwLength;
				state.lastTag = element.dwTag;
				state.elemIt++;
			}
			continue;
		}

		if (state.seq == &m_RootElement) {
			state.seq->dwLength = pItem->dwLength;
			break;
		}

		for (;;) {
			SequenceItem *pItem = *state.itemIt;
			pItem->dwHeaderSize = CalculateItemHeaderSize(pItem->dwLength);
			state.seq->dwLength += pItem->dwHeaderSize + pItem->dwLength;

			state.itemIt++;
			if (state.itemIt == state.seq->items.end()) {
				Element *pSeq = state.seq;
				state = stateStack.top();
				stateStack.pop();
				pSeq->dwHeaderSize = CalculateElementHeaderSize(pSeq->dwTag, state.lastTag, pSeq->dwLength);
				SequenceItem *pItem = *state.itemIt;
				pItem->dwLength += pSeq->dwHeaderSize + pSeq->dwLength;
				//state.seq->dwLength += pSeq->dwHeaderSize + pSeq->dwLength;
				state.elemIt++;
				break;
			} else {
				if (pItem->elements.size() != 0 && pItem->fragmented) {
					pItem->dwLength = 0;
					state.elemIt = (*state.itemIt)->elements.begin();
					state.lastTag = 0;
					break;
				}
			}
		}
	}

	// allocate new buffer
	dwDataBufferLength = m_RootElement.dwLength;
	pDataBuffer = new BYTE[m_RootElement.dwLength];
	if (pDataBuffer == NULL)
		return false;
	buf.pCursor = pDataBuffer;
	buf.dwRemaining = dwDataBufferLength;

	state.seq = &m_RootElement;
	state.itemIt = m_RootElement.items.begin();
	state.elemIt = (*state.itemIt)->elements.begin();
	state.lastTag = 0;

	// write elements
	for (;;) {
		if (state.elemIt != (*state.itemIt)->elements.end()) {
			Element &element = state.elemIt->second;
			element.fragmented = false;
			if (!PutElementHeader(buf, state.lastTag, &element)) {
				delete pDataBuffer;
				return false;
			}
			state.lastTag = element.dwTag;
			if (element.type == etSequence && element.items.size() != 0) {
				element.pData = buf.pCursor;
				element.dataAllocated = false;
				stateStack.push(state);
				state.seq = &element;
				state.itemIt = state.seq->items.begin();
				state.elemIt = (*state.itemIt)->elements.begin();
				state.lastTag = 0;
				(*state.itemIt)->fragmented = false;
				if (!PutItemHeader(buf, *state.itemIt)) {
					delete pDataBuffer;
					return false;
				}
				(*state.itemIt)->pData = buf.pCursor;
			} else {
				BYTE *pData = buf.pCursor;
				if (element.dwLength != 0) {
					if (!PutBytes(buf, element.dwLength, element.pData)) {
						delete pDataBuffer;
						return false;
					}
				}
				if (element.dataAllocated)
					delete[] element.pData;
				element.pData = pData;
				element.dataAllocated = false;
				state.elemIt++;
			}
			continue;
		}

		if (state.seq == &m_RootElement) {
			break;
		}

		for (;;) {
			state.itemIt++;
			if (state.itemIt == state.seq->items.end()) {
				state = stateStack.top();
				stateStack.pop();
				state.elemIt++;
				break;
			} else {
				(*state.itemIt)->fragmented = false;
				if (!PutItemHeader(buf, *state.itemIt)) {
					delete pDataBuffer;
					return false;
				}
				(*state.itemIt)->pData = buf.pCursor;
				if ((*state.itemIt)->elements.size() != 0) {
					state.elemIt = (*state.itemIt)->elements.begin();
					state.lastTag = 0;
					break;
				}
			}
		}
	}

	delete m_RootElement.pData;

	m_RootElement.fragmented = false;
	m_RootElement.items[0]->fragmented = false;
	m_RootElement.pData = pDataBuffer;
	m_RootElement.items[0]->pData = pDataBuffer;
	m_RootElement.dwLength = dwDataBufferLength;
	m_RootElement.items[0]->dwLength = dwDataBufferLength;

	return true;
}


DWORD CTaggedDataBuffer::CalculateElementHeaderSize(DWORD dwTag, DWORD dwLastTag, DWORD dwDataLength)
{
	DWORD dwHeaderSize;

	dwHeaderSize = 1;

	dwTag -= dwLastTag;
	if (dwTag <= 0xFFUL)
		dwHeaderSize += 1;
	else if (dwTag <= 0xFFFFUL)
		dwHeaderSize += 2;
	else if (dwTag <= 0xFFFFFFUL)
		dwHeaderSize += 3;
	else
		dwHeaderSize += 4;

	if (dwDataLength < 1 || dwDataLength > 28) {
		if (dwDataLength <= 0xFFUL)
			dwHeaderSize += 1;
		else if (dwDataLength <= 0xFFFFUL)
			dwHeaderSize += 2;
		else if (dwDataLength <= 0xFFFFFFUL)
			dwHeaderSize += 3;
		else
			dwHeaderSize += 4;
	}

	return dwHeaderSize;
}


DWORD CTaggedDataBuffer::CalculateItemHeaderSize(DWORD dwItemLength)
{
	if (dwItemLength <= 252)
		return 1;
	else if (dwItemLength <= 0xFFFFUL)
		return 3;
	else if (dwItemLength <= 0xFFFFFFUL)
		return 4;
	else
		return 5;
}


bool CTaggedDataBuffer::PutElementHeader(BufferState &bufState, DWORD dwLastTag, Element *pElement)
{
	DWORD dwElementLength = pElement->dwLength;
	BYTE hdr[9];
	DWORD dwHdrLen;
	DWORD dwHdrTag;

	hdr[0] = 0;
	dwHdrLen = 1;

	if (pElement->type == etSequence)
		hdr[0] |= 0x80;

	dwHdrTag = pElement->dwTag - dwLastTag;
	if (dwHdrTag <= 0xFFUL) {
		hdr[dwHdrLen] = (BYTE)dwHdrTag;
		dwHdrLen += 1;
	} else if (dwHdrTag <= 0xFFFFUL) {
		hdr[0] |= 0x20;
		hdr[dwHdrLen] = (BYTE)dwHdrTag;
		hdr[dwHdrLen+1] = (BYTE)(dwHdrTag >> 8);
		dwHdrLen += 2;
	} else if (dwHdrTag <= 0xFFFFFFUL) {
		hdr[0] |= 0x40;
		hdr[dwHdrLen] = (BYTE)dwHdrTag;
		hdr[dwHdrLen+1] = (BYTE)(dwHdrTag >> 8);
		hdr[dwHdrLen+2] = (BYTE)(dwHdrTag >> 16);
		dwHdrLen += 3;
	} else {
		hdr[0] |= 0x60;
		hdr[dwHdrLen] = (BYTE)dwHdrTag;
		hdr[dwHdrLen+1] = (BYTE)(dwHdrTag >> 8);
		hdr[dwHdrLen+2] = (BYTE)(dwHdrTag >> 16);
		hdr[dwHdrLen+3] = (BYTE)(dwHdrTag >> 24);
		dwHdrLen += 4;
	}

	if (dwElementLength >= 1 && dwElementLength <= 28) {
		hdr[0] |= dwElementLength - 1;
	} else if (dwElementLength <= 0xFFUL) {
		hdr[0] |= 28;
		hdr[dwHdrLen] = (BYTE)dwElementLength;
		dwHdrLen += 1;
	} else if (dwElementLength <= 0xFFFFUL) {
		hdr[0] |= 29;
		hdr[dwHdrLen] = (BYTE)dwElementLength;
		hdr[dwHdrLen+1] = (BYTE)(dwElementLength >> 8);
		dwHdrLen += 2;
	} else if (dwElementLength <= 0xFFFFFFUL) {
		hdr[0] |= 30;
		hdr[dwHdrLen] = (BYTE)dwElementLength;
		hdr[dwHdrLen+1] = (BYTE)(dwElementLength >> 8);
		hdr[dwHdrLen+2] = (BYTE)(dwElementLength >> 16);
		dwHdrLen += 3;
	} else {
		hdr[0] |= 31;
		hdr[dwHdrLen] = (BYTE)dwElementLength;
		hdr[dwHdrLen+1] = (BYTE)(dwElementLength >> 8);
		hdr[dwHdrLen+2] = (BYTE)(dwElementLength >> 16);
		hdr[dwHdrLen+3] = (BYTE)(dwElementLength >> 24);
		dwHdrLen += 4;
	}

	return PutBytes(bufState, dwHdrLen, hdr);
}


bool CTaggedDataBuffer::PutItemHeader(BufferState &bufState, SequenceItem *pItem)
{
	DWORD dwItemLength = pItem->dwLength;
	BYTE hdr[5];
	DWORD dwHdrLen;

	if (dwItemLength <= 252) {
		hdr[0] = (BYTE)dwItemLength;
		dwHdrLen = 1;
	} else if (dwItemLength <= 0xFFFFUL) {
		hdr[0] = 253;
		hdr[1] = (BYTE)dwItemLength;
		hdr[2] = (BYTE)(dwItemLength >> 8);
		dwHdrLen = 3;
	} else if (dwItemLength <= 0xFFFFFFUL) {
		hdr[0] = 254;
		hdr[1] = (BYTE)dwItemLength;
		hdr[2] = (BYTE)(dwItemLength >> 8);
		hdr[3] = (BYTE)(dwItemLength >> 16);
		dwHdrLen = 4;
	} else {
		hdr[0] = 255;
		hdr[1] = (BYTE)dwItemLength;
		hdr[2] = (BYTE)(dwItemLength >> 8);
		hdr[3] = (BYTE)(dwItemLength >> 16);
		hdr[3] = (BYTE)(dwItemLength >> 24);
		dwHdrLen = 5;
	}

	return PutBytes(bufState, dwHdrLen, hdr);
}

