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


int checkCondition(const wchar_t *pszCondition, bool bComplex)
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
	yystate.scanBufferLength = (bUnicode ? len / sizeof(wchar_t) : len) + 3;
	yystate.scanBuffer = new wchar_t[yystate.scanBufferLength];
	if (!yystate.scanBuffer)
	{
		fclose(in);
		fwprintf(stderr, L"fatal error: Cannot allocate buffer for source file: '%s'\n", pszFileName);
		return NULL;
	}
	memset(yystate.scanBuffer, 0, yystate.scanBufferLength * sizeof(wchar_t));

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
		for (size_t pos = 0, idx = 0; pos < len / sizeof(wchar_t); pos++, idx++)
		{
			wchar_t c = _mbbtombc(fgetc(in));
			if (c == '\r')
			{
				if (idx < len / sizeof(wchar_t) - 1)
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
