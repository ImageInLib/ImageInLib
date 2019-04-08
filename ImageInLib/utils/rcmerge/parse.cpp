#include "headers.h"
#include "grammar.h"
#ifndef fileno
#define fileno _fileno
#endif
#ifndef read
#define read _read
#endif
#ifndef isatty
#define isatty _isatty
#endif
#include "parse.h"
#include "rcfile.h"
#include <mbstring.h>


#define REPLACEMENT_CHAR 0xFFFD

// Unicode surrogate pairs constants
#define UNI_SUR_HIGH_START	0xD800
#define UNI_SUR_HIGH_END	0xDBFF
#define UNI_SUR_LOW_START	0xDC00
#define UNI_SUR_LOW_END		0xDFFF
#define UNI_SUR_MASK		0x3FF
#define UNI_SUR_SHIFT		10


extern FILE* yyin;

extern int yyparse(void *yystate);
extern void yyrestart(FILE*input_file);
extern struct yy_buffer_state * yy_scan_buffer(wchar_t *base, unsigned int size);
extern void yy_switch_to_buffer(yy_buffer_state * new_buffer);


void parseError(yy_parse_state *yystate, const char *pszErrorMsg)
{
	fflush(stderr);
	fwprintf(stderr, L"%s(%d,%d) : error : %S\n", yystate->fileName, yystate->curLine, yystate->curColumn, pszErrorMsg);
}


int  checkCondition(const wchar_t *pszCondition, bool bComplex)
{
	if (bComplex) {
		if (wcsstr(pszCondition, L"!defined(_WIN32)") != NULL)
			return 0;
		else if (wcsstr(pszCondition, L"defined(_WIN32)") != NULL)
			return 1;
		if (wcsstr(pszCondition, L"!defined(APSTUDIO_INVOKED)") != NULL)
			return 0;
		else if (wcsstr(pszCondition, L"defined(APSTUDIO_INVOKED)") != NULL)
			return 1;
	} else {
		if (wcscmp(pszCondition, L"_WIN32") == 0)
			return 1;
		if (wcscmp(pszCondition, L"APSTUDIO_INVOKED") == 0)
			return 1;
	}
	return -1;
}


RCFile_C *parseRCFile(const wchar_t *pszFileName, bool bRequireUnicode)
{
	yy_buffer_state * yybufstate;
	yy_parse_state yystate;
	int res;

	yyin = NULL;

	// open input .RC file
	FILE *in = NULL;
	if (_wfopen_s(&in, pszFileName, L"rb") != 0)
	{
		fwprintf(stderr, L"fatal error: Cannot open source file: '%s'\n", pszFileName);
		return NULL;
	}

	// detect file character set
	bool bUnicode = false;
	char c1, c2;
	c1 = fgetc(in);
	c2 = fgetc(in);
	fseek(in, 0, SEEK_SET);
	size_t len = _filelength(_fileno(in));
	if (len >= 2 && c1 == '\xff' && c2 == '\xfe') 
	{
		// Unicode with BOM
		bUnicode = true;
	}
	else if (len >= 2 && c1 != 0 && c2 == 0) 
	{
		// Unicode without BOM
		bUnicode = true;
	}
	else
	{
		// 8-bit
		bUnicode = false;
	}

	// check encoding
	if (bRequireUnicode && !bUnicode)
	{
		fwprintf(stderr, L"fatal error: Invalid encoding of source file (only Unicode is accepted): '%s'\n", pszFileName);
		return NULL;
	}

	// allocate parser buffer
	yystate.scanBufferLength = (bUnicode ? len / sizeof(wchar_t) : len) + 2;
	yystate.scanBuffer = new wchar_t[yystate.scanBufferLength + 3];
	if (!yystate.scanBuffer)
	{
		fclose(in);
		fwprintf(stderr, L"fatal error: Cannot allocate buffer for source file: '%s'\n", pszFileName);
		return NULL;
	}
	memset(yystate.scanBuffer, 0, (yystate.scanBufferLength + 3) * sizeof(wchar_t));

	// read .RC file
	if (bUnicode)
	{
		for (size_t pos = 0, idx = 0; pos < len / sizeof(wchar_t); pos++)
		{
			wchar_t c = fgetwc(in);
			if (c == L'\xfeff')
			{
				yystate.scanBufferLength--;
			}
			else if (c == '\r')
			{
				if (idx < len / sizeof(wchar_t) - 1)
				{
					c = fgetwc(in);
					if (c != '\n')
						yystate.scanBuffer[idx++] = '\r';
					else
						yystate.scanBufferLength--;
					pos++;
				}
				yystate.scanBuffer[idx++] = c;
			}
			else
			{
				yystate.scanBuffer[idx++] = c;
			}
		}
	}
	else
	{
		for (size_t pos = 0, idx = 0; pos < len; pos++, idx++)
		{
			wchar_t c = _mbbtombc(fgetc(in));
			if (c == '\r')
			{
				if (idx < len - 1)
				{
					c = _mbbtombc(fgetc(in));
					if (c != '\n')
						yystate.scanBuffer[idx++] = '\r';
					else
						yystate.scanBufferLength--;
					pos++;
				}
			}
			yystate.scanBuffer[idx] = c;
		}
	}

	// close .RC file
	fclose(in);
	in = NULL;

	// set up scanner buffer state
	yybufstate  = yy_scan_buffer(yystate.scanBuffer, yystate.scanBufferLength);

	// initialize parser state
    yystate.curTokenStart = 0;
    yystate.curTokenEnd = 0;
    yystate.curLine = 1;
    yystate.curColumn = 1;
    yystate.trigSkip = 0;
	yystate.fileName = pszFileName;
	yystate.rc = new RCFile_C;
	yystate.conditionStack.push(-1);
	yystate.conditionState = -1;
	yystate.stringtable = NULL;

	// parse .RC file
	res = yyparse(&yystate);

	if (res != 0) {
		delete yystate.rc;
		return NULL;
    }
	return yystate.rc;
}


std::wstring quoteString(const wchar_t *str)
{
	std::wstring s;
	s = L"\"";
	if (*str == '\'' && *(str+1) == '-') {
		s += L"-";
		str++;
	}
	while (*str) {
		if (*str == '\"')
			s += L"\"\"";
		else
			s += *str;
		str++;
	}
	s += L"\"";
	return s;
}


wchar_t linebuf[65536];

RCFile_C *parseTXTFile(const wchar_t *pszFileName, bool bRequireUnicode)
{
	RCFile_C *rc = NULL;

	// open input .TXT file
	FILE *in = NULL;
	if (_wfopen_s(&in, pszFileName, L"rb") != 0)
	{
		fwprintf(stderr, L"fatal error: Cannot open source file: '%s'\n", pszFileName);
		return NULL;
	}

	// detect file character set
	bool bUnicode = false;
	bool bUTF8 = false;
	char c1, c2, c3;
	c1 = fgetc(in);
	c2 = fgetc(in);
	c3 = fgetc(in);
	fseek(in, 0, SEEK_SET);
	size_t len = _filelength(_fileno(in));
	if (len >= 3 && c1 == '\xef' && c2 == '\xbb' && c3 == '\xbf') 
	{
		// UTF-8 with BOM
		bUTF8 = true;
		len--;
	}
	else if (len >= 2 && c1 == '\xff' && c2 == '\xfe') 
	{
		// Unicode with BOM
		bUnicode = true;
		len--;
	}
	else if (len >= 2 && c1 != 0 && c2 == 0) 
	{
		// Unicode without BOM
		bUnicode = true;
	}
	else
	{
		// 8-bit
		bUnicode = false;
	}

	// check encoding
	if (bRequireUnicode && !bUnicode)
	{
		fwprintf(stderr, L"fatal error: Invalid encoding of source file (only Unicode is accepted): '%s'\n", pszFileName);
		return NULL;
	}

	// allocate parser buffer
	size_t scanBufferLength;
	wchar_t * scanBuffer;
	scanBufferLength = (bUnicode ? len / sizeof(wchar_t) : len) + 3;
	scanBuffer = new wchar_t[scanBufferLength];
	if (!scanBuffer)
	{
		fclose(in);
		fwprintf(stderr, L"fatal error: Cannot allocate buffer for source file: '%s'\n", pszFileName);
		return NULL;
	}
	memset(scanBuffer, 0, scanBufferLength * sizeof(wchar_t));

	// read .TXT file
	if (bUnicode)
	{
		for (size_t pos = 0, idx = 0; pos < len / sizeof(wchar_t); pos++)
		{
			wchar_t c = fgetwc(in);
			if (c == L'\xfeff')
			{
				scanBufferLength--;
			}
			else if (c == '\r')
			{
				if (idx < len / sizeof(wchar_t) - 1)
				{
					c = fgetwc(in);
					if (c != '\n')
						scanBuffer[idx++] = '\r';
					else
						scanBufferLength--;
					pos++;
				}
				scanBuffer[idx++] = c;
			}
			else
			{
				scanBuffer[idx++] = c;
			}
		}
	}
	else if (bUTF8)
	{
		for (size_t pos = 0, idx = 0; pos < len; pos++)
		{
			// UTF-8 to UCS-2
			unsigned char c;
			wchar_t w;
			c = fgetc(in);
			if (c < 0x80)
			{
				if (c == '\r')
				{
					if (idx < len - 1)
					{
						c = fgetc(in);
						if (c != '\n')
							scanBuffer[idx++] = '\r';
						else
							scanBufferLength--;
						pos++;
					}
				}
				scanBuffer[idx++] = c;
			}
			else
			{
				if (c <= 0xBF)
				{
					scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR; // invalid character sequence
				}
				else if (c <= 0xDF)
				{ 
					// 110XXXXX 10XXXXXX
					if (len == 0)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					w = ((c & 0x1F) << 6);
					c = fgetc(in);
					--len;
					if (c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					w += (c & 0x3F);
					if (w <= 0x7F) {
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					scanBuffer[idx++] = w;
				} 
				else if (c <= 0xEF)
				{ 
					// 1110XXXX 10XXXXXX 10XXXXXX
					if (len == 0)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					w = ((c & 0x0F) << 12);
					c = fgetc(in);
					--len;
					if (len == 0 || c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					w += ((c & 0x3F) << 6);
					c = fgetc(in);
					--len;
					if (c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					w += (c & 0x3F);
					if (w <= 0x7FF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					scanBuffer[idx++] = w;
				}
				else if (c <= 0xF7)
				{ 
					// 11110XXX 10XXXXXX 10XXXXXX 10XXXXXX
					int utf32;
					if (len == 0)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					utf32 = ((c & 0x07) << 18);
					c = fgetc(in);
					--len;
					if (len == 0 || c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					utf32 += ((c & 0x3F) << 12);
					c = fgetc(in);
					--len;
					if (len == 0 || c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					utf32 += ((c & 0x3F) << 6);
					c = fgetc(in);
					--len;
					if (c < 0x80 || c > 0xBF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					utf32 += (c & 0x3F);
					if (utf32 <= 0xFFFF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					if (utf32 > 0x10FFFF)
					{
						scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR;
						continue; // invalid character sequence
					}
					// 0xFFFF - 0x10FFFF
					utf32 -= 0x10000L;
					scanBuffer[idx++] = (wchar_t)((utf32 >> UNI_SUR_SHIFT) + UNI_SUR_HIGH_START);
					scanBuffer[idx++] = (wchar_t)((utf32 & UNI_SUR_MASK) + UNI_SUR_LOW_START);
				}
				else
				{
					scanBuffer[idx++] = (wchar_t)REPLACEMENT_CHAR; // invalid character sequence
				}
			}
		}
	}
	else
	{
		for (size_t pos = 0, idx = 0; pos < len; pos++, idx++)
		{
			wchar_t c = _mbbtombc(fgetc(in));
			if (c == '\r')
			{
				if (idx < len - 1)
				{
					c = _mbbtombc(fgetc(in));
					if (c != '\n')
						scanBuffer[idx++] = '\r';
					else
						scanBufferLength--;
					pos++;
				}
			}
			scanBuffer[idx] = c;
		}
	}

	// close .TXT file
	fclose(in);
	in = NULL;

	Dialog_C *pDialog = NULL;
	Menu_C *pMenu = NULL;
	rc = new RCFile_C();
	StringTable_C *pStringTable = new StringTable_C(rc, L"");
	rc->AddResource(pStringTable);
	std::map<std::wstring, MenuPopup_C*> popups;

	// load .TXT file
	int line = 1;
	size_t scanBufferPos = 0;
	while (scanBufferPos != scanBufferLength)
	{
		std::wstring type, id, orig, trans;

		linebuf[0] = 0;
		linebuf[_countof(linebuf)-1] = 0;
		for (size_t idx = 0; scanBufferPos < scanBufferLength && idx < _countof(linebuf)-1; idx++, scanBufferPos++)
		{
			linebuf[idx] = scanBuffer[scanBufferPos];
			if(linebuf[idx] == '\n')
			{
				linebuf[idx+1] = 0;
				scanBufferPos++;
				break;
			}
		}
		if (*linebuf)
		{
			int len = wcslen(linebuf);
			while (len && (linebuf[len-1] == '\r' || linebuf[len-1] == '\n'))
				linebuf[--len] = 0;
		}

		wchar_t *p = linebuf;
		while (*p && *p != '\t')
			type += *p++;
		if (!*p)
		{
			fprintf(stdout, "FATAL ERROR: invalid file format at line %d\n", line);
			goto fail;
		}
		p++;
		while (*p && *p != '\t')
			id += *p++;
		if (!*p)
		{
			fprintf(stdout, "FATAL ERROR: invalid file format at line %d\n", line);
			goto fail;
		}
		p++;
		while (*p && *p != '\t')
			orig += *p++;
		if (!*p)
		{
			fprintf(stdout, "FATAL ERROR: invalid file format at line %d\n", line);
			goto fail;
		}
		p++;
		trans = p;

		if (type == L"STRING")
		{
			pDialog = NULL;
			pMenu = NULL;
			popups.clear();
			String_C *pString = new String_C(id.c_str(), quoteString(trans.c_str()).c_str(), 0, 0);
			pString->SetTranslationHint(quoteString(orig.c_str()).c_str());
			pStringTable->AddString(pString); 
		}
		else if (type == L"DIALOG")
		{
			pMenu = NULL;
			popups.clear();
			pDialog = new Dialog_C(rc, id.c_str(), trans.size() ? quoteString(trans.c_str()).c_str() : L"" , 0, 0);
			pDialog->SetTranslationHint(orig.size() ? quoteString(orig.c_str()).c_str() : L"");
			rc->AddResource(pDialog);
		}
		else if (type == L"MENU")
		{
			pDialog = NULL;
			pMenu = new Menu_C(rc, id.c_str(), L"");
			popups.clear();
			popups[id.c_str()] = NULL;
			rc->AddResource(pMenu);
		}
		else if (pDialog != NULL)
		{
			const wchar_t *ctlid = id.c_str();
			if (wcsncmp(ctlid, pDialog->GetName(), wcslen(pDialog->GetName())) != 0)
			{
				fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
				goto fail;
			}
			ctlid += wcslen(pDialog->GetName());
			if (*ctlid != ':')
			{
				fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
				goto fail;
			}
			ctlid++;
			Control_C *pControl = new Control_C(type.c_str(), quoteString(trans.c_str()).c_str(), ctlid, L"", 0, 0);
			pControl->SetTranslationHint(quoteString(orig.c_str()).c_str());
			pDialog->AddControl(pControl);
		}
		else if (pMenu != NULL && type == L"POPUP")
		{
			const wchar_t *p = wcsrchr(id.c_str(), L':');
			if (p == NULL)
			{
				fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
				goto fail;
			}
			else
			{
				std::map<std::wstring, MenuPopup_C*>::iterator it;
				it = popups.find(std::wstring(id.c_str(), p - id.c_str()));
				if (it == popups.end()) {
					fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
					goto fail;
				}
				MenuPopup_C *pPopup = new MenuPopup_C(it->second, quoteString(trans.c_str()).c_str(), *(p+1) == '#' ? L"" : p+1, L"", 0, 0);
				pPopup->SetTranslationHint(quoteString(orig.c_str()).c_str());
				if (it->second == NULL)
					pMenu->AddItem(pPopup);
				else
					it->second->AddItem(pPopup);
				popups[id.c_str()] = pPopup;
			}
		}
		else if (pMenu != NULL && type == L"MENUITEM")
		{
			const wchar_t *p = wcsrchr(id.c_str(), L':');
			if (p == NULL)
			{
				fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
				goto fail;
			}
			else
			{
				std::map<std::wstring, MenuPopup_C*>::iterator it;
				it = popups.find(std::wstring(id.c_str(), p - id.c_str()));
				if (it == popups.end())
				{
					fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
					goto fail;
				}
				MenuItem_C *pItem = new MenuItem_C(quoteString(trans.c_str()).c_str(), p+1, L"", 0, 0);
				pItem->SetTranslationHint(quoteString(orig.c_str()).c_str());
				if (it->second == NULL)
					pMenu->AddItem(pItem);
				else
					it->second->AddItem(pItem);
			}
		}
		else
		{
			fwprintf(stdout, L"FATAL ERROR: invalid file format at line %d\n", line);
			goto fail;
		}

		line++;
	}

	// close .RC file
	return rc;

fail:
	if (rc)
		delete rc;
	return NULL;
}
