#include "headers.h"
#include "parse.h"
#include "rcfile.h"
#include "typeinfo.h"

static char *usage_text = 
	"Usage: " __APPNAME__ " [options] input.rc <build> [<compared_idr>]\n"
	"\n"
	"options:\n"
	"  -U | --unicode  accept only files with Unicode charset\n"
#ifdef _DEBUG
	"  --yytrace   traces internal parser states\n"
#endif
	"  --help        display this help and exit\n"
	"  --version     output version information and exit\n"
	"  --show-match  show matching items\n"
	"build:\n"
	"  TOMOCON_WORKSTATION | TOMOCON_LITE | TOMOCON_VIEWER | TOMOCON_DEBUG\n"
	"compared_idr:\n"
	"  IDs of menus compared to IDR_HOTKEY\n"
	"\n"
	"prints two column comparison. left column consists of ids in IDR_HOTKEY\n"
	"\n";

static void usage() 
{
	printf(usage_text);
}


static const char *version_text = 
	__APPNAME__ " version %s\n"
	"Copyright (C) 2011 TatraMed s.r.o.\n";

static void version() 
{
	printf(version_text, __APP_VERSION__);
}

bool getCommands( std::set<std::wstring> & commands, Menu_C * pMenu, const wchar_t filter)
{
	if ( pMenu == NULL)
		return false;

	stack<MenuItem_C *> items;

	for ( int i = 0; i < pMenu->GetItemsCount(); ++i)
		items.push( pMenu->GetItem( i) );

	for ( ; items.size() > 0; )
	{
		MenuItem_C * pItem = items.top();
		items.pop();

		if ( MenuPopup_C * pPopup = dynamic_cast<MenuPopup_C *>( pItem) )
		{
			for ( int i = 0; i < pPopup->GetItemsCount(); i++)
				items.push( pPopup->GetItem( i));
		}
		else if ( pItem != NULL && wcscmp( pItem->GetID(), L"SEPARATOR") != 0)
		{
			if ( filter != 0)
			{
				wchar_t * p = (wchar_t*) wcsstr( pItem->GetText(), L"~");
				
				if ( p == NULL)
				{
					commands.insert( pItem->GetID());
				}
				else
					for ( p++; *p != 0; ++p)
				{
					if ( *p == filter)
					{
						commands.insert( pItem->GetID());
						break;
					}
				}				
			}
			else
			{
				commands.insert( pItem->GetID());
			}
		}
	}

	return true;
}

bool compareSets( std::set<std::wstring> & setA, std::set<std::wstring> & setB, bool bHideMatch)
{
	bool bRet = true;

	std::set<std::wstring>::iterator ii = setA.begin();
	std::set<std::wstring>::iterator jj = setB.begin();

	for ( ; ii != setA.end() && jj != setB.end(); )
	{
		if ( *ii < *jj)
		{
			bRet = false;
			wprintf( L"%s\n", (*ii).c_str());
			++ii;
		}
		else if ( *ii > *jj)
		{
			bRet = false;
			printf( "                                    %s\n", (*jj).c_str());
			++jj;
		}
		else  // *ii == *jj
		{
			if ( ! bHideMatch)
				wprintf( L"%-35s %s\n", (*ii).c_str(), (*jj).c_str());
			++jj;
			++ii;
		}
	}

	for ( ; jj != setB.end(); ++jj)
	{
		bRet = false;
		wprintf( L"                                    %s\n", (*jj).c_str());
	}
	for ( ; ii != setA.end(); ++ii)
	{
		bRet = false;
		wprintf( L"%s\n", (*ii).c_str());
	}

	return bRet;
}

int wmain(int argc, wchar_t **argv)
{	
	bool bRequireUnicode = false;
	bool bHideMatch = true;
	const wchar_t *rcname;
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
		}
		else if ( ! wcscmp( argv[i], L"--show-match") )
		{
			bHideMatch = false;
#ifdef _DEBUG
		}
		else if (!wcscmp(argv[i], L"--yytrace"))
		{
			extern int yydebug;
			yydebug = 1;
#endif
		}
		else if ( (argc - i) < 3)
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

	wchar_t * builds[] = { L"TOMOCON_WORKSTATION", L"TOMOCON_RELEASE", L"TOMOCON_VIEWER", L"TOMOCON_LITE", L"TOMOCON_DEBUG", NULL};
	wchar_t flags[] = L"WWVLD";
	int nBuild = 0;

	++i;
	for ( ; builds[nBuild] != NULL; ++nBuild)
	{
		if( ! wcscmp( builds[nBuild], argv[i]))
			break;
	}

	std::set<std::wstring> idrs;

	for ( ++i; i < argc; ++i)
		idrs.insert( argv[i]);

	if ( builds[ nBuild] == NULL || idrs.empty() )
	{
		usage();
		exit(1);
	}

	std::set<std::wstring> buildset;
	buildset.insert( builds[ nBuild]);
	if ( ! wcscmp( builds[ nBuild], L"TOMOCON_WORKSTATION") )
		buildset.insert( L"TOMOCON_RELEASE");
	else if ( ! wcscmp( builds[ nBuild], L"TOMOCON_RELEASE") )
		buildset.insert( L"TOMOCON_WORKSTATION");
	else if ( ! wcscmp( builds[ nBuild], L"TOMOCON_DEBUG") )
		buildset.insert( L"TOMOCON_WORKSTATION");

	std::set<std::wstring> hotkeyIDs;
	std::set<std::wstring> menusIDs;

	fwprintf(stdout, L"Parsing resources...\n");

	fwprintf(stdout, L"%s\n", rcname);

	rc = parseRCFile(rcname, bRequireUnicode);
	if (rc == NULL)
		goto fail;

	for ( int resIdx = 0; resIdx < rc->GetResourceCount(); resIdx++) 
		if ( Menu_C * pMenu = dynamic_cast<Menu_C*>( rc->GetResource( resIdx)) )
	{
		std::wstring name = pMenu->GetName();
		std::wstring def = L"";

		int i = name.find(L"$(");
		int j = name.rfind(L")");
		if ( i != string::npos && j != string::npos)
		{
			def = name.substr( i + 2, j - i - 2);
			name.resize( i);
		}

		if ( def.empty() == false && buildset.find( def) == buildset.end() )
			continue;	// not wished build

		if ( name == L"IDR_HOTKEY")
		{ 
			getCommands( hotkeyIDs, pMenu, flags[nBuild]);
		}
		else
		{
			std::set<std::wstring>::iterator ii = idrs.find( name);
			if ( ii != idrs.end() )
			{
				getCommands( menusIDs, pMenu, 0);
			}
		}
	}

	compareSets( hotkeyIDs, menusIDs, bHideMatch);

	if (rc)
		delete rc;

    return 0;

fail:
	if (rc)
		delete rc;
	return 1;
}
