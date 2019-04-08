#include "headers.h"
#include "parse.h"
#include "merge.h"
#include "rcfile.h"

bool g_verbose = false;
std::wstring g_mark;

static wchar_t *usage_text = 
	L"Usage: " __APPNAME__ L" [options] input.rc translation output.rc\n"
	L"\n"
	L"options:\n"
	L"  --mark <text>   mark missing or different entries with <text> mark\n"
#ifdef _DEBUG
	L"  --yytrace   traces internal parser states\n"
#endif
	L"  --verbose   report also missing translations\n"
	L"  --help      display this help and exit\n"
	L"  --version   output version information and exit\n"
	L"  --rc		   load translation in .RC format\n";

static void usage() 
{
	wprintf(usage_text);
}


static const wchar_t *version_text = 
	__APPNAME__ L" version %s\n"
	L"Copyright (C) 2004 TatraMed s.r.o.\n";

static void version() 
{
	wprintf(version_text, __APP_VERSION__);
}


int wmain(int argc, wchar_t **argv)
{	
	bool bRequireUnicode = false;
	const wchar_t *rcname, *trans, *mergerc;
	RCFile_C *rc = NULL;
	RCFile_C *xrc = NULL;
	bool loadRCTranslation = false;
	int i;

	if (argc <= 1)
	{
		usage();
		exit(1);
	}
	for(i = 1; i < argc; i++)
	{
		if (!wcscmp(argv[i], L"-U") || !wcscmp(argv[i], L"--unicode"))
		{
			bRequireUnicode = true;
		}
		if (!wcscmp(argv[i], L"--help")) 
		{
			usage();
			exit(1);
		}
		else if (!wcscmp(argv[i], L"--version"))
		{
			version();
			exit(1);
		}
		else if (!wcscmp(argv[i], L"--mark"))
		{
			if (argc - i < 4) {
				usage();
				exit(1);
			}
			g_mark = argv[++i];
		}
		else if (!wcscmp(argv[i], L"--verbose") || !wcscmp(argv[i], L"-v"))
		{
			g_verbose = true;
		}
		else if (!wcscmp(argv[i], L"--rc"))
		{
			loadRCTranslation = true;
#ifdef _DEBUG
		}
		else if (!wcscmp(argv[i], L"--yytrace"))
		{
			extern int yydebug;
			yydebug = 1;
#endif
		}
		else if (argc - i != 3)
		{
			usage();
			exit(1);
		}
		else
		{
			rcname = argv[i];
			trans = argv[i+1];
			mergerc = argv[i+2];
			break;
		}
	}

	fwprintf(stdout, L"Parsing resources...\n");
	fwprintf(stdout, L"%s\n", rcname);
	rc = parseRCFile(rcname, bRequireUnicode);
	if (rc == NULL)
		goto fail;

	fwprintf(stdout, L"Parsing translated resources...\n");
	if (loadRCTranslation)
	{
		fwprintf(stdout, L"%s\n", trans);
		xrc = parseRCFile(trans, bRequireUnicode);
		if (xrc == NULL)
			goto fail;
	}
	else 
	{
		fwprintf(stdout, L"%s\n", trans);
		xrc = parseTXTFile(trans, bRequireUnicode);
		if (xrc == NULL)
			goto fail;
	}

	fwprintf(stdout, L"Merging resources...\n");
	mergeTexts(rc, xrc, rcname, mergerc);

	if (rc)
		delete rc;
	if (xrc)
		delete xrc;

    return 0;

fail:
	if (rc)
		delete rc;
	if (xrc)
		delete xrc;
	return 1;
}
