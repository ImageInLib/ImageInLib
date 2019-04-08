// TaggedDataBuffer.h: interface for the CTaggedDataBuffer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TAGGEDDATABUFFER_H__22109CA3_99FF_4ECE_BAED_6BF04B313B90__INCLUDED_)
#define AFX_TAGGEDDATABUFFER_H__22109CA3_99FF_4ECE_BAED_6BF04B313B90__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class _HTDBELEMENT {
public:
	virtual ~_HTDBELEMENT() {}
};
class _HTDBITEM {
public:
	virtual ~_HTDBITEM() {}
};

typedef _HTDBELEMENT *HTDBELEMENT;
typedef _HTDBITEM *HTDBITEM;

class CTaggedDataBuffer  
{
public:
	CTaggedDataBuffer();
	virtual ~CTaggedDataBuffer();

	bool Parse(const BYTE *pBuffer, DWORD dwBufferSize);
	bool Write(BYTE *pBuffer, DWORD *pdwBufferSize);

	void Clear();

	enum ElementType {
		etValue,
		etSequence
	};

	HTDBITEM GetItem(HTDBELEMENT hElement, DWORD dwIndex);
	HTDBITEM GetNextItem(HTDBITEM hItem);
	DWORD GetItemsCount(HTDBELEMENT hElement);

	HTDBELEMENT GetFirstElement(HTDBITEM hItem);
	HTDBELEMENT GetNextElement(HTDBELEMENT hElement);

	HTDBELEMENT FindElement(HTDBITEM hItem, DWORD dwTag);

	bool GetItemInfo(HTDBITEM hItem, DWORD* pdwTag, DWORD* pdwIndex, DWORD* pdwLength, BYTE** ppData);
	bool GetElementInfo(HTDBELEMENT hElement, ElementType* pType, DWORD* pdwTag, DWORD* pdwLength, BYTE** ppData);

	bool SetElementValue(HTDBITEM hItem, ElementType nType, DWORD dwTag, DWORD dwLength, const BYTE* pData);
	HTDBELEMENT CreateSequence(HTDBITEM hItem, DWORD dwTag);
	HTDBITEM CreateSequenceItem(HTDBELEMENT hSequence);

// attributes
protected:
	class SequenceItem;
	class Element;

	struct BufferState {
		BYTE *pCursor;
		DWORD dwRemaining;
	};

	typedef vector<SequenceItem*> SequenceItems;

	class Element : public _HTDBELEMENT {
	public:
		Element() : dataAllocated(false) {}
		~Element() { 
			for_each(items.begin(), items.end(), DeleteSequenceItem); 
			if (dataAllocated) 
				delete pData; 
		}
	public:
		ElementType type;
		DWORD dwTag;
		DWORD dwHeaderSize;
		DWORD dwLength;
		BYTE *pData;
		SequenceItems items;
		Element *pParent;
		DWORD itemIdx;
		bool dataAllocated;
		bool fragmented;
	protected:
		static void DeleteSequenceItem(class SequenceItem *pItem) { delete pItem; }
	};

	typedef map<DWORD, Element> ElementMap;

	class SequenceItem : public _HTDBITEM {
	public:
		DWORD dwHeaderSize;
		DWORD dwLength;
		BYTE *pData;
		DWORD itemIdx;
		ElementMap elements;
		Element *pParent;
		bool fragmented;
	};

	struct ParserState {
		BufferState	buf;
		ElementMap *elements;
		DWORD lastTag;
		Element *seq;
		BufferState	seqBuf;
	};

	struct DefragState {
		Element *seq;
		ElementMap::iterator elemIt;
		SequenceItems::iterator itemIt;
		DWORD lastTag;
	};

	Element m_RootElement;

// implementation
protected:

	void Init();

	bool ParseSequence(Element *pSequence);

	bool ParseElement(BufferState &bufState, ElementType &elementType,
		DWORD &headerSize, DWORD &elementTag, DWORD &elementLength, BYTE* &pElementData);

	bool ParseSequenceItem(BufferState &bufState, DWORD &headerSize, DWORD &itemLength, BYTE*& pItemData);

	bool DefragmentBuffer();

	DWORD CalculateElementHeaderSize(DWORD dwTag, DWORD dwLastTag, DWORD dwDataLength);
	DWORD CalculateItemHeaderSize(DWORD dwItemLength);

	bool PutElementHeader(BufferState &bufState, DWORD dwLastTag, Element *pElement);
	bool PutItemHeader(BufferState &bufState, SequenceItem *pItem);

	inline DWORD GetAvailable(BufferState &bufState) {
		return bufState.dwRemaining;
	}

	inline bool IsAvailable(BufferState &bufState, DWORD dwBytes) {
		return bufState.dwRemaining >= dwBytes;
	}

	inline bool GetByte(BufferState &bufState, BYTE &b) {
		if (bufState.dwRemaining != 0) {
			b = *bufState.pCursor++;
			bufState.dwRemaining--;
			return true;
		}
		return false;
	}

	inline bool GetBytes(BufferState &bufState, DWORD dwBytes, BYTE* &pBuffer) {
		if (bufState.dwRemaining >= dwBytes) {
			pBuffer = bufState.pCursor;
			bufState.pCursor += dwBytes;
			bufState.dwRemaining -= dwBytes;
			return true;
		}
		return false;
	}

	inline bool PutBytes(BufferState &bufState, DWORD dwBytes, BYTE* pBuffer) {
		if (bufState.dwRemaining >= dwBytes) {
			memmove(bufState.pCursor, pBuffer, dwBytes);
			bufState.pCursor += dwBytes;
			bufState.dwRemaining -= dwBytes;
			return true;
		}
		return false;
	}

};

#endif // !defined(AFX_TAGGEDDATABUFFER_H__22109CA3_99FF_4ECE_BAED_6BF04B313B90__INCLUDED_)
