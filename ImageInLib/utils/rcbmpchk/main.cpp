#include "headers.h"
#include "parse.h"
#include "check.h"
#include "rcfile.h"

static char *usage_text = 
	"Usage: " __APPNAME__ " [options] input.rc\n"
	"\n"
	"options:\n"
	"  -U | --unicode  accept only files with Unicode charset\n"
#ifdef _DEBUG
	"  --yytrace       traces internal parser states\n"
#endif
	"  --help          display this help and exit\n"
	"  --basedir DIRECTORY\n"
	"                  opens bitmap files relative to this directory\n"
	"  --version       output version information and exit\n";

static void usage() 
{
	printf(usage_text);
}


static const char *version_text = 
	__APPNAME__ " version %s\n"
	"Copyright (C) 2011 TatraMed Software s.r.o.\n";

static void version() 
{
	printf(version_text, __APP_VERSION__);
}


int wmain(int argc, wchar_t **argv)
{	
	bool bCheckBitmaps = true;
	bool bCheckCursors = true;
	bool bRequireUnicode = false;
	int res = 0;
	const wchar_t *rcname = L"";
	const wchar_t *basedir = L"";
	RCFile_C *rc = NULL;
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
		else if (!wcscmp(argv[i], L"--help"))
		{
			usage();
			exit(1);
		}
		else if (!wcscmp(argv[i], L"--version"))
		{
			version();
			exit(1);
		}
#ifdef _DEBUG
		else if (!wcscmp(argv[i], L"--yytrace"))
		{
			extern int yydebug;
			yydebug = 1;
		}
#endif
		else if (!wcscmp(argv[i], L"--basedir"))
		{
			if (argc - i < 1) 
			{
				usage();
				exit(1);
			}
			basedir = argv[++i];
		}
		else if (!wcscmp(argv[i], L"--no-bitmaps"))
		{
			bCheckBitmaps = false;
		}
		else if (!wcscmp(argv[i], L"--no-cursors"))
		{
			bCheckCursors = false;
		}
		else if (argc - i != 1)
		{
			usage();
			exit(1);
		}
		else
		{
			rcname = argv[i];
			break;
		}
	}

	fwprintf(stdout, L"Parsing resources...\n");

	fwprintf(stdout, L"%s\n", rcname);
	rc = parseRCFile(rcname, bRequireUnicode);
	if (rc == NULL)
		goto fail;

	if (bCheckBitmaps)
	{
		fwprintf(stdout, L"Checking bitmaps...\n");
		if (checkBitmaps(rc, basedir))
			res = 1;
	}

	if (bCheckCursors)
	{
		fwprintf(stdout, L"Checking cursors...\n");
		if (checkCursors(rc, basedir))
			res = 1;
	}

fail:
	if (rc)
		delete rc;

    return res;
}
