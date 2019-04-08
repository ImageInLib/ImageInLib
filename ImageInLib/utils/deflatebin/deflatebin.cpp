#include <stdio.h>
#include <tchar.h>
#include <io.h>
#include <string.h>
#include <stdlib.h>
#include "zlib/zlib.h"
#include <Windows.h>

_TCHAR use [] = _T("bin2deflate <inputfile> <output.h> <variablename>\n");

#define LINE_WIDTH 16

int _tmain(int argc, _TCHAR* argv[])
{
	unsigned char *inputBuffer;
	unsigned long inputLength;
	unsigned char *outputBuffer;
	unsigned long outputLength;

	if ( argc < 4)
	{
		_tprintf(_T("%s: wrong arguments !\n"), *argv, use);
		_tprintf(use);
		exit(1);
	}

	FILE  *in;
	if (_tfopen_s(&in, argv[1], _T("rb")) != 0)
	{
		_tprintf(_T("error: failed to open %s\n"), argv[1]);
		exit(1);
	}

	inputLength = _filelength(_fileno(in));
	if (inputLength > 256*1024*1024)
	{
		_tprintf(_T("error: input file is too long\n"));
		exit(1);
	}

	inputBuffer = (unsigned char*)malloc(inputLength);
	if (!inputBuffer)
	{
		_tprintf(_T("error: could not allocate buffer for input file\n"));
		exit(1);
	}

	if (fread(inputBuffer, 1, inputLength, in) != inputLength)
	{
		_tprintf(_T("error: could not read input file\n"));
		exit(1);
	}

	outputLength = inputLength + inputLength/10 + 1024;
	outputBuffer = (unsigned char*)malloc(outputLength);
	if (!outputBuffer)
	{
		_tprintf(_T("error: could not allocate buffer for output file\n"));
		exit(1);
	}
	if (compress2(outputBuffer, &outputLength, inputBuffer, inputLength, Z_BEST_COMPRESSION) != Z_OK)
	{
		_tprintf(_T("error: could not compress input file\n"));
		exit(1);
	}

	FILE *out;
	if (_tfopen_s(&out, argv[2], _T("wt")) != 0)
	{
		_tprintf(_T("error: failed to create %s\n"), argv[2]);
		exit(1);
	}

	if (_ftprintf(out, _T("#pragma once\n\nstatic unsigned char %s [] = \n{\n"), argv[3]) < 0)
	{
		_tprintf(_T("error: failed write to %s !\n"), argv[0], argv[1]);
		exit(1);
	}
	if (_ftprintf(out, _T("\t0x%02x,0x%02x,0x%02x,0x%02x,\n"), 
		(inputLength & 0xff), ((inputLength >> 8) & 0xff), ((inputLength >> 16) & 0xff), ((inputLength >> 24) & 0xff)) < 0)
	{
		_tprintf(_T("error: failed write to %s !\n"), argv[0], argv[1]);
		exit(1);
	}
	if (_ftprintf(out, _T("\t0x%02x,0x%02x,0x%02x,0x%02x,\n"), 
		(outputLength & 0xff), ((outputLength >> 8) & 0xff), ((outputLength >> 16) & 0xff), ((outputLength >> 24) & 0xff)) < 0)
	{
		_tprintf(_T("error: failed write to %s !\n"), argv[0], argv[1]);
		exit(1);
	}
	unsigned char *p = outputBuffer;
	while((p - outputBuffer) < (ptrdiff_t)outputLength)
	{
		if (_ftprintf(out, _T("\t")) < 0)
		{
			_tprintf(_T("error: failed write to %s !\n"), argv[2]);
			exit(1);
		}
		for (unsigned long cnt = min(LINE_WIDTH, outputBuffer + outputLength - p); cnt; --cnt, ++p)
		{
			if (_ftprintf(out, _T("0x%02x%s"), *p, p != outputBuffer + outputLength - 1 ? "," : "") < 0)
			{
				_tprintf(_T("error: failed write to %s !\n"), argv[2]);
				exit(1);
			}
		}
		if (_ftprintf(out, _T("\n")) < 0)
		{
			_tprintf(_T("error: failed write to %s !\n"), argv[2]);
			exit(1);
		}
	}
	if (_ftprintf(out, _T("};\n")) < 0)
	{
		_tprintf(_T("error: failed write to %s !\n"), argv[2]);
		exit(1);
	}

	return 0;
}

