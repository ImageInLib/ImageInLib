#include <stdio.h>
#include <tchar.h>

_TCHAR use [] = _T("%s [-c] <output.h> <variablename> <fname1> [<fname2> [<fname3>]...]\n\
\tall input files are merged into one array\n\n");

#define LINE_WIDTH 16

int _tmain(int argc, _TCHAR* argv[])
{
	bool charType = false;

	if ( argc > 1 && _tcscmp(argv[1], _T("-c")) == 0)
	{
		charType = true;
		--argc;
		for (int i = 1; i < argc; ++i)
			argv[i] = argv[i+1];
	}

	if ( argc < 4)
	{
		_tprintf( _T("%s: wrong number of arguments!\n"), argv[0]);
		_tprintf( use, *argv);
		return -1;
	}
	FILE  * F;
	if (_tfopen_s(&F, argv[1], _T("wt")) != 0)
	{
		_tprintf( _T("%s: failed to create %s\n"), argv[0], argv[1]);
		return -1;
	}

	if ( _ftprintf( F, _T("#pragma once\n\nstatic const %schar %s [] = \n{\n"), charType ? _T("") : _T("unsigned "), argv[2]) < 0)
	{
		_tprintf( _T("%s: failed to write %s !\n"), argv[0], argv[1]);
		fclose( F);
		return -1;
	}

	int ret = 0;
	size_t cnt = 0;

	unsigned char buffer[LINE_WIDTH];
	for ( int i = 3; i < argc && ret == 0; ++i)
	{
		FILE * G = NULL;
		if (_tfopen_s(&G, argv[i], _T("rb")) != 0) {
			ret = -1;
			_tprintf( _T("%s: failed to open %s\n"), argv[0], argv[i]);
			continue;
		}
		for ( ; feof( G) == 0;)
		{
			size_t sz = fread( buffer, sizeof( unsigned char), LINE_WIDTH, G);
			if ( sz < LINE_WIDTH && ferror( G) != 0)
			{
				_tprintf( _T("%s: failed to read %s ! \n"), argv[0], argv[i]);
				ret = -1;
				break;
			}

			if ( _ftprintf( F, _T("\t")) < 0)
			{
				_tprintf( _T("%s: failed to write %s !\n"), argv[0], argv[1]);
				ret = -1;
				break;
			}
			for ( size_t j = 0; j < sz; ++j)
			{
				if ( _ftprintf( F, _T("0x%02x,"), buffer[j]) < 0)
				{
					_tprintf( _T("%s: failed to write %s !\n"), argv[0], argv[1]);
					ret = -1;
					break;
				}
			}
			if ( _ftprintf( F, _T("\n")) < 0)
			{
				_tprintf( _T("%s: failed to write %s !\n"), argv[0], argv[1]);
				ret = -1;
				break;
			}

			cnt += sz;
		}
		fclose( G);
	}

	if (ret == 0)
	{
		if ( _ftprintf( F, _T("\t0x00\n};\n\nstatic const unsigned int %s_len = %u;\n\n"), argv[2], cnt) < 0)
		{
			_tprintf( _T("%s: failed to write %s !\n"), argv[0], argv[1]);
			ret = -1;
		}
	}

	_tprintf( _T("%s: created '%s' with array '%s' of %u bytes\n"), argv[0], argv[1], argv[2], cnt);

	fclose( F);
	return ret;
}

