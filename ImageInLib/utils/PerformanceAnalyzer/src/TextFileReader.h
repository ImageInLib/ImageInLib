#pragma once
class CTextFileReader
{
public:
	CTextFileReader(LPCTSTR lpszFileName);
	~CTextFileReader();

	bool OpenFile();
	bool CloseFile();
	int  ReadLine(CString& strLine);

protected:
	int ReadToBuffer(CString& strBuffer);

protected:
	CString m_strFileName;
	FILE*   m_file;
	CString m_strBuffer;
	int     m_iLastPosition;
};

