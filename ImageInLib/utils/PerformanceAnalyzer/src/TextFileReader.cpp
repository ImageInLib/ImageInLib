#include "stdafx.h"
#include "TextFileReader.h"


CTextFileReader::CTextFileReader(LPCTSTR lpszFileName) :
	m_strFileName(lpszFileName),
	m_file(nullptr),
	m_iLastPosition(0)
{
}


CTextFileReader::~CTextFileReader()
{
	if (m_file)
		fclose(m_file);
}

bool CTextFileReader::OpenFile()
{
	return (_tfopen_s(&m_file, m_strFileName, _T("rt")) == 0);
}

bool CTextFileReader::CloseFile()
{
	int nResult = fclose(m_file);
	m_file = nullptr;
	return (nResult == 0);
}

int CTextFileReader::ReadLine(CString& strLine)
{
	// try to find end of line in the buffer string
	int iPosition = m_strBuffer.Find(TCHAR('\n'), m_iLastPosition);
	// not found, try to read from file
	if (iPosition == -1) {
		// delete data that was already read
		m_strBuffer.Delete(0, m_iLastPosition);
		m_iLastPosition = 0;
		// read new data
		int nResult = ReadToBuffer(m_strBuffer);
		if (nResult != 0)
			return nResult;
		if (m_strBuffer.IsEmpty())
			return -1;	// end of file reached
		// try to find end of line in the buffer string again. if not found, set position to the end.
		iPosition = m_strBuffer.Find(TCHAR('\n'), m_iLastPosition);
		if (iPosition == -1)
			iPosition = m_strBuffer.GetLength();
	}
	// copy the line to strLine
	_ASSERTE(iPosition >= m_iLastPosition);
	strLine = m_strBuffer.Mid(m_iLastPosition, iPosition - m_iLastPosition);
	m_iLastPosition = iPosition + 1;
	return 0;
}

int CTextFileReader::ReadToBuffer(CString& strBuffer)
{
	const UINT nBufferLength = 65536;
	CHAR buffer[nBufferLength];

	size_t nResult;
	nResult = fread_s(buffer, nBufferLength - 1, sizeof(CHAR), nBufferLength - 1, m_file);
	if (nResult == 0 && ferror(m_file) != 0)
		return errno;
	buffer[nResult] = 0;

	strBuffer += CString(buffer);
	return 0;
}
