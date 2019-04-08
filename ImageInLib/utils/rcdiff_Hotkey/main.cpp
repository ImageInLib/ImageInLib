#include "headers.h"
#include "parse.h"
#include "compare.h"
#include "rcfile.h"

static char *usage_text = 
	"Usage: " __APPNAME__ " [options] file1.rc file2.rc\n"
	"\n"
	"\t- compares IDR_HOTKEY menu items and their ~ flags\n"
	"\t- checks strings for IDR_HOTKEY menu in file2.rc\n"
	"\n"
	"options:\n"
	"  -U | --unicode  accept only files with Unicode charset\n"
#ifdef _DEBUG
	"  --yytrace       traces internal parser states\n"
#endif
	"  --help          display this help and exit\n"
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
	bool bRequireUnicode = false;
	const wchar_t *rc1name = L"", *rc2name = L"";
	RCFile_C *rc1 = NULL, *rc2 = NULL;
	int i;

	for(i = 1; i < argc; i++) {
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
#ifdef _DEBUG
		}
		else if (!wcscmp(argv[i], L"--yytrace"))
		{
			extern int yydebug;
			yydebug = 1;
#endif
		} 
		else if (argc - i != 2)
		{
			usage();
			exit(1);
		}
		else
		{
			rc1name = argv[i];
			rc2name = argv[i+1];
			break;
		}
	}

	fwprintf(stdout, L"Parsing resources...\n");

	fwprintf(stdout, L"%s\n", rc1name);
	rc1 = parseRCFile(rc1name, bRequireUnicode);
	if (rc1 == NULL)
		goto fail;

	fwprintf(stdout, L"%s\n", rc2name);
	rc2 = parseRCFile(rc2name, bRequireUnicode);
	if (rc2 == NULL)
		goto fail;

	fwprintf(stdout, L"Comparing resources...\n");

	if (compareFiles(rc1, rc2) == 0) {
		fwprintf(stdout, L"no differences found.\n");
		fwprintf(stdout, L"\n");
	} else {
		fwprintf(stdout, L"\n");
		goto fail;
	}

	if (rc1)
		delete rc1;
	if (rc2)
		delete rc2;

    return 0;

fail:
	if (rc1)
		delete rc1;
	if (rc2)
		delete rc2;
	return 1;
}
