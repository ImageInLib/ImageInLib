#pragma once

#include <Windows.h>
#include <list>
#include <algorithm>
#include <limits.h>

class DiBitsEncoder
{
private:
	DiBitsEncoder(const DiBitsEncoder&);            // private undefined copy constructor	
	DiBitsEncoder& operator=(const DiBitsEncoder&); // private undefined copy assignment operator

	static const size_t BlockSize = 65536;

public:

	DiBitsEncoder()
		: m_nWorkBufferPos(0)
		, m_pOutputBlock(new BYTE[DiBitsEncoder::BlockSize])
		, m_OutputBlockPos(0)
	{
		memset(m_WorkBuffer, 0, sizeof(m_WorkBuffer));
	}

	~DiBitsEncoder()
	{
		std::for_each(m_OutputBlocks.begin(), m_OutputBlocks.end(), [](BYTE *buf) { delete[] buf; });
		delete[] m_pOutputBlock;
	}

	inline void Encode(const BYTE *pBuffer, size_t nNumberOfBytes)
	{
		if (!pBuffer)
			return;

		while (nNumberOfBytes) 
		{
			int nLengthCode;
			size_t nBlockSize;
			int nEncodedLength;
			int nDiBitOrder[4];

			// find best block size
			nEncodedLength = INT_MAX;
			nBlockSize = 256;
			for (int nTestLengthCode = 15; nTestLengthCode >= 0; nTestLengthCode--)
			{
				size_t nTestBlockSize;
				int nDiBitCounts[4];
				int nTestDiBitOrder[4];
				int nTestEncodedLength;

				// get block size:
				// 0x0*-0x9* = 1-10 bytes
				// 0xA* = 16 bytes
				// 0xB* = 32 bytes
				// 0xC* = 64 bytes
				// 0xD* = 128 bytes
				// 0xE* = 192 bytes
				// 0xF* = 256 bytes
				if (nTestLengthCode == 0xf)
        			nTestBlockSize = 0x100;
				else if (nTestLengthCode == 0xe)
        			nTestBlockSize = 0xc0;
				else if (nTestLengthCode >= 0xa)
        			nTestBlockSize = 0x10 << (nTestLengthCode - 0xa);
				else
        			nTestBlockSize = nTestLengthCode + 1;

				if (nTestBlockSize > nNumberOfBytes)
					continue;

				// get di-bit counts
				nDiBitCounts[3] = nDiBitCounts[2] = nDiBitCounts[1] = nDiBitCounts[0] = 0;
				for (size_t i = 0; i < nTestBlockSize; i++)
				{
					int b = pBuffer[i];
					nDiBitCounts[b & 3]++;
					nDiBitCounts[(b >> 2) & 3]++;
					nDiBitCounts[(b >> 4) & 3]++;
					nDiBitCounts[(b >> 6) & 3]++;
				}

				// sort di-bits
				nTestDiBitOrder[0] = 0;
				nTestDiBitOrder[1] = 1;
				nTestDiBitOrder[2] = 2;
				nTestDiBitOrder[3] = 3;
				if (nDiBitCounts[nTestDiBitOrder[0]] < nDiBitCounts[nTestDiBitOrder[1]])
					std::swap(nTestDiBitOrder[0], nTestDiBitOrder[1]);
				if (nDiBitCounts[nTestDiBitOrder[0]] < nDiBitCounts[nTestDiBitOrder[2]])
					std::swap(nTestDiBitOrder[0], nTestDiBitOrder[2]);
				if (nDiBitCounts[nTestDiBitOrder[0]] < nDiBitCounts[nTestDiBitOrder[3]])
					std::swap(nTestDiBitOrder[0], nTestDiBitOrder[3]);
				if (nDiBitCounts[nTestDiBitOrder[1]] < nDiBitCounts[nTestDiBitOrder[2]])
					std::swap(nTestDiBitOrder[1], nTestDiBitOrder[2]);
				if (nDiBitCounts[nTestDiBitOrder[1]] < nDiBitCounts[nTestDiBitOrder[3]])
					std::swap(nTestDiBitOrder[1], nTestDiBitOrder[3]);

				// di-bit type 2 must be less or equal to di-bit type 3
				if (nTestDiBitOrder[2] > nTestDiBitOrder[3])
					std::swap(nTestDiBitOrder[2], nTestDiBitOrder[3]);

				// calculate encoded length
				nTestEncodedLength = nDiBitCounts[nTestDiBitOrder[0]] +
									 2 * nDiBitCounts[nTestDiBitOrder[1]] +
									 3 * nDiBitCounts[nTestDiBitOrder[2]] +
									 3 * nDiBitCounts[nTestDiBitOrder[3]];
				nTestEncodedLength = ((nTestEncodedLength + 7) >> 3) + 1;

				if ((double)nTestEncodedLength / nTestBlockSize < (double)nEncodedLength/nBlockSize)
				{
					nLengthCode = nTestLengthCode;
					nBlockSize = nTestBlockSize;
					nEncodedLength = nTestEncodedLength;
					nDiBitOrder[0] = nTestDiBitOrder[0];
					nDiBitOrder[1] = nTestDiBitOrder[1];
					nDiBitOrder[2] = nTestDiBitOrder[2];
					nDiBitOrder[3] = nTestDiBitOrder[3];
				}
//				else
//				{
//					// compression ratio for this block size is worse than previous, we have a winner
//					break; 
//				}
			}

			// encode di-bits block
			int bTypeAndSize = (nLengthCode << 4) | (nDiBitOrder[0] << 2) | nDiBitOrder[1];
			m_WorkBuffer[(m_nWorkBufferPos >> 3)] |= (BYTE)(bTypeAndSize << (m_nWorkBufferPos & 7));
			m_WorkBuffer[(m_nWorkBufferPos >> 3) + 1] |= (BYTE)((bTypeAndSize << (m_nWorkBufferPos & 7)) >> 8);
			m_nWorkBufferPos += 8;
			for (size_t i = 0; i < nBlockSize; i++)
			{
				int b = pBuffer[i];
				int bitbuf = 0;
				int bitcnt = 0;
				for (int j = 0; j < 4; j++)
				{
					int dibit = b & 3;
					b >>= 2;
					if (dibit == nDiBitOrder[0])
					{
						bitbuf = (bitbuf << 1) | 1;	// 1
						bitcnt++;
					}
					else if (dibit == nDiBitOrder[1])
					{
						bitbuf = (bitbuf << 2);		// 00
						bitcnt += 2;
					}
					else if (dibit == nDiBitOrder[2])
					{
						bitbuf = (bitbuf << 3) | 2;	// 010
						bitcnt += 3;
					}
					else
					{
						bitbuf = (bitbuf << 3) | 6;	// 011
						bitcnt += 3;
					}
				}
				m_WorkBuffer[(m_nWorkBufferPos >> 3)] |= (BYTE)(bitbuf << (m_nWorkBufferPos & 7));
				m_WorkBuffer[(m_nWorkBufferPos >> 3) + 1] |= (BYTE)((bitbuf << (m_nWorkBufferPos & 7)) >> 8);
				m_WorkBuffer[(m_nWorkBufferPos >> 3) + 2] |= (BYTE)((bitbuf << (m_nWorkBufferPos & 7)) >> 16);
				m_nWorkBufferPos += bitcnt;
			}

			for (size_t i=0; i < (m_nWorkBufferPos >> 3); i++)
			{
				if (m_OutputBlockPos == DiBitsEncoder::BlockSize)
				{
					m_OutputBlocks.push_back(m_pOutputBlock);
					m_pOutputBlock = new unsigned char[DiBitsEncoder::BlockSize];
					m_OutputBlockPos = 0;
				}
				m_pOutputBlock[m_OutputBlockPos++] = m_WorkBuffer[i];
			}
			m_WorkBuffer[0] = m_WorkBuffer[m_nWorkBufferPos >> 3];
			memset(m_WorkBuffer + 1, 0, sizeof(m_WorkBuffer) - 1);
			m_nWorkBufferPos &= 7;

			nNumberOfBytes -= nBlockSize;
			pBuffer += nBlockSize;
		}
	}

	inline size_t GetOutputSize()
	{
		size_t result = m_OutputBlockPos + DiBitsEncoder::BlockSize * m_OutputBlocks.size() + ((m_nWorkBufferPos + 7) >> 3) + 1;
		return result;
	}

	inline void WriteOutput(BYTE *pDestBuffer)
	{
		if (pDestBuffer)
		{
			std::for_each(m_OutputBlocks.begin(), m_OutputBlocks.end(), [&pDestBuffer](BYTE *pSrcBuffer) { 
				memcpy(pDestBuffer, pSrcBuffer, DiBitsEncoder::BlockSize);
				pDestBuffer += DiBitsEncoder::BlockSize;
			});
			m_WorkBuffer[(m_nWorkBufferPos >> 3)] |= (BYTE)(0xff << (m_nWorkBufferPos & 7));
			m_WorkBuffer[(m_nWorkBufferPos >> 3) + 1] |= (BYTE)((0xff << (m_nWorkBufferPos & 7)) >> 8);
			m_nWorkBufferPos += 8;
			for (size_t i=0; i < ((m_nWorkBufferPos + 7) >> 3); i++)
			{
				if (m_OutputBlockPos == DiBitsEncoder::BlockSize)
				{
					m_OutputBlocks.push_back(m_pOutputBlock);
					m_pOutputBlock = new unsigned char[DiBitsEncoder::BlockSize];
					m_OutputBlockPos = 0;
				}
				m_pOutputBlock[m_OutputBlockPos++] = m_WorkBuffer[i];
			}
			m_nWorkBufferPos = 0;
			if (m_OutputBlockPos > 0)
			{
				memcpy(pDestBuffer, m_pOutputBlock, m_OutputBlockPos);
				pDestBuffer += m_OutputBlockPos;
			}
		}
	}

private:
	BYTE   m_WorkBuffer[512];           // work buffer of di-bits encoder
	size_t m_nWorkBufferPos;            // number of bits in work buffer
	BYTE   *m_pOutputBlock;				// current output block
	size_t  m_OutputBlockPos;			// current position in output block
	std::list<BYTE*> m_OutputBlocks;	// list of full output blocks
};

