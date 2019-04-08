#include "headers.h"
#include "parse.h"
#include "extract.h"
#include "rcfile.h"

static char *usage_text = 
	"Usage: " __APPNAME__ " [options] input.rc output.txt\n"
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
	"Copyright (C) 20011 TatraMed s.r.o.\n";

static void version() 
{
	printf(version_text, __APP_VERSION__);
}


int wmain(int argc, wchar_t **argv)
{	
	bool bRequireUnicode = false;
	const wchar_t *rcname = L"", *txtname = L"";
	RCFile_C *rc = NULL;
	int i;

	if (argc <= 1) {
		usage();
		exit(1);
	}
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
			rcname = argv[i];
			txtname = argv[i+1];
			break;
		}
	}

	fwprintf(stdout, L"Parsing resources...\n");

	fwprintf(stdout, L"%s\n", rcname);
	rc = parseRCFile(rcname, bRequireUnicode);
	if (rc == NULL)
		goto fail;

	fwprintf(stdout, L"Extracting texts...\n");
	extractTexts(rc, txtname);

	if (rc)
		delete rc;

    return 0;

fail:
	if (rc)
		delete rc;
	return 1;
}
