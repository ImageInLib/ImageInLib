#pragma once

#include <Windows.h>
#include <list>
#include <algorithm>

class RLEEncoder
{
private:
	RLEEncoder(const RLEEncoder&);            // private undefined copy constructor	
	RLEEncoder& operator=(const RLEEncoder&); // private undefined copy assignment operator

	static const size_t BlockSize = 65536;

public:

	RLEEncoder()
		: m_RepeatBufferPos(1)
		, m_PrevLiteral(-1)
		, m_RepeatCount(0)
		, m_pOutputBlock(new BYTE[RLEEncoder::BlockSize])
		, m_OutputBlockPos(0)
	{
		m_RepeatBuffer[0] = 0;
	}

	~RLEEncoder()
	{
		std::for_each(m_OutputBlocks.begin(), m_OutputBlocks.end(), [](BYTE *buf) { delete[] buf; });
		delete[] m_pOutputBlock;
	}

	inline void Encode(const BYTE *pBuffer, size_t nNumberOfBytes)
	{
		if (!pBuffer)
			return;

		while (nNumberOfBytes--) 
		{
			int b = (int)*pBuffer++;

			// if the current byte equals the last byte read
			// (which is initialized with the "impossible" value -1),
			// just increase the repeat counter
			if (m_PrevLiteral == (int)b) 
			{
				m_RepeatCount++;
			}
			else
			{
				// byte is different from last byte read.
				// flush replicate run if necessary
				switch (m_RepeatCount)
				{
				case 0:
					// empty buffer
					break;

				case 1:
					// literal run of one byte
					m_RepeatBuffer[m_RepeatBufferPos++] = (BYTE)m_PrevLiteral;
					break;

				case 2:
					// literal run of two bytes
					m_RepeatBuffer[m_RepeatBufferPos++] = (BYTE)m_PrevLiteral;
					m_RepeatBuffer[m_RepeatBufferPos++] = (BYTE)m_PrevLiteral;
					break;

				default:
					// more than two bytes in repeat buffer. Convert to replicate run
					if (m_RepeatBufferPos > 1)
					{
						// there is a literal run in the buffer that must be flushed
						// before the replicate run.  Flush literal run now.
						m_RepeatBuffer[0] = (BYTE)(m_RepeatBufferPos-2);
						WriteEncodedBytes(m_RepeatBufferPos);
					}

					// this is the byte value for the repeat run
					m_RepeatBuffer[1] = (BYTE)m_PrevLiteral;

					// write as many repeat runs as necessary
					for (; m_RepeatCount>0; m_RepeatCount-=128)
					{
						if (m_RepeatCount > 128) 
							m_RepeatBuffer[0] = 129;
						else m_RepeatBuffer[0] = 
							(BYTE)(257 - m_RepeatCount);
						WriteEncodedBytes(2);
					}

					// now the buffer is guaranteed to be empty
					m_RepeatBuffer[0] = 0;
					m_RepeatBufferPos = 1;
					break;
				}

				// if we have 128 or more bytes in the literal run, flush repeat buffer
				if (m_RepeatBufferPos > 129)
				{
					m_RepeatBuffer[0] = 127;
					WriteEncodedBytes(129);
					m_RepeatBufferPos -= 128;
					if (m_RepeatBufferPos > 1)
						m_RepeatBuffer[1] = m_RepeatBuffer[129];
					if (m_RepeatBufferPos > 2)
						m_RepeatBuffer[2] = m_RepeatBuffer[130];
				}

				// current byte is stored in m_PrevLiteral, m_RepeatCount is 1.
				m_PrevLiteral = b;
				m_RepeatCount = 1;
			}
		}
	}

	inline size_t GetOutputSize()
	{
		Finish();
		size_t result = m_OutputBlockPos + RLEEncoder::BlockSize * m_OutputBlocks.size();
		return result;
	}

	inline void WriteOutput(BYTE *pDestBuffer)
	{
		if (pDestBuffer)
		{
			Finish();
			std::for_each(m_OutputBlocks.begin(), m_OutputBlocks.end(), [&pDestBuffer](BYTE *pSrcBuffer) { 
				memcpy(pDestBuffer, pSrcBuffer, RLEEncoder::BlockSize);
				pDestBuffer += RLEEncoder::BlockSize;
			});
			if (m_OutputBlockPos > 0)
				memcpy(pDestBuffer, m_pOutputBlock, m_OutputBlockPos);
		}
	}

protected:
	inline void WriteEncodedBytes(size_t numberOfBytes)
	{
		for (size_t i=0; i < numberOfBytes; i++)
		{
			if (m_OutputBlockPos == RLEEncoder::BlockSize)
			{
				m_OutputBlocks.push_back(m_pOutputBlock);
				m_pOutputBlock = new unsigned char[RLEEncoder::BlockSize];
				m_OutputBlockPos = 0;
			}
			m_pOutputBlock[m_OutputBlockPos++] = m_RepeatBuffer[i];
		}
	}

	inline void Finish()
	{
		// if there are max 1 bytes in the repeat counter, convert to literal run
		if (m_RepeatCount < 2)
		{
			for (; m_RepeatCount > 0; --m_RepeatCount) 
				m_RepeatBuffer[m_RepeatBufferPos++] = (BYTE)m_PrevLiteral;
		}

		// if we have 128 or more bytes in the literal run, flush buffer
		if (m_RepeatBufferPos > 129)
		{
			m_RepeatBuffer[0] = 127;
			WriteEncodedBytes(129);
			m_RepeatBufferPos -= 128;
			if (m_RepeatBufferPos > 1)
				m_RepeatBuffer[1] = m_RepeatBuffer[129];
			if (m_RepeatBufferPos > 2)
				m_RepeatBuffer[2] = m_RepeatBuffer[130];
		}

		// if there is still a literal run in the buffer, flush literal run
		if (m_RepeatBufferPos > 1)
		{
			m_RepeatBuffer[0] = (BYTE)(m_RepeatBufferPos-2);
			WriteEncodedBytes(m_RepeatBufferPos);
		}

		// if there is a remaining repeat run, flush this one as well
		if (m_RepeatCount >= 2)
		{
			m_RepeatBuffer[1] = (BYTE)m_PrevLiteral;
			// write as many repeat runs as necessary
			for (; m_RepeatCount>0; m_RepeatCount-=128)
			{
				if (m_RepeatCount > 128)
					m_RepeatBuffer[0] = 0x81;
				else 
					m_RepeatBuffer[0] = (BYTE)(257 - m_RepeatCount);
				WriteEncodedBytes(2);
			}
		}

		// now the buffer is guaranteed to be empty, re-initialize
		m_RepeatBuffer[0] = 0;
		m_PrevLiteral = -1;
		m_RepeatCount = 0;
		m_RepeatBufferPos = 1;
	}

private:
	BYTE	m_RepeatBuffer[132];		// work buffer of RLE encoder
	size_t	m_RepeatBufferPos;			// current position in work buffer
	int		m_PrevLiteral;				// previous byte added
	int		m_RepeatCount;				// bytes repeat count

	BYTE   *m_pOutputBlock;				// current output block
	size_t  m_OutputBlockPos;			// current position in output block
	std::list<BYTE*> m_OutputBlocks;	// list of full output blocks
};

