#include "headers.h"
#include "compare.h"
#include "rcfile.h"
#include "typeinfo.h"

typedef std::map<std::wstring, Resource_C*> resources_map;
typedef std::map<std::wstring, std::wstring> string_map;
typedef std::list<std::string> string_list;

int compareMenu(Menu_C *menu1, Menu_C *menu2, string_map & strMap)
{
	MenuPopup_C *pPopup1, *pPopup2;
	stack<int> popupStack;
	int itemIdx;
	int nRet = 0;

	if (menu1->GetItemsCount() != menu2->GetItemsCount()) {
		wprintf(L"differences in menu: %s\n", menu1->GetName());
		return 1;
	}

	// compare menu
	pPopup1 = NULL;
	pPopup2 = NULL;
	itemIdx = 0;
	for (;;) {
		MenuItem_C *pItem1, *pItem2;
		if (pPopup1)
			pItem1 = pPopup1->GetItem(itemIdx);
		else
			pItem1 = menu1->GetItem(itemIdx);
		if (pPopup2)
			pItem2 = pPopup2->GetItem(itemIdx);
		else
			pItem2 = menu2->GetItem(itemIdx);
		itemIdx++;
		if (pItem1 == NULL) {
			if (popupStack.size() == 0)
				break;
			itemIdx = popupStack.top();
			popupStack.pop();
			pPopup1 = pPopup1->GetOwner();
			pPopup2 = pPopup2->GetOwner();
			continue;
		}
		if (pItem1 == NULL || pItem2 == NULL || wcscmp(pItem1->GetID(), pItem2->GetID()) != 0) {
			wprintf(L"differences in menu: %s\n", menu1->GetName());
			return 1;
		}
		if ((dynamic_cast<MenuPopup_C*>(pItem1) != NULL) != (dynamic_cast<MenuPopup_C*>(pItem2) != NULL)) {
			wprintf(L"differences in menu: %s\n", menu1->GetName());
			return 1;
		}
		if (dynamic_cast<MenuPopup_C*>(pItem1)) {
			popupStack.push(itemIdx);
			pPopup1 = dynamic_cast<MenuPopup_C*>(pItem1);
			pPopup2 = dynamic_cast<MenuPopup_C*>(pItem2);
			if (pPopup1->GetItemsCount() != pPopup1->GetItemsCount()) {
				wprintf(L"differences in menu: %s\n", menu1->GetName());
				return 1;
			}
			itemIdx = 0;
		}
		else
		{	
			// compare post ~ flags as order independent
			const wchar_t *p1 = pItem1->GetText();
			const wchar_t *p2 = pItem2->GetText();
			
			if ( p1 != NULL)
				for ( ; *p1 != 0 && *p1 != '~'; )
					p1++;
			
			if ( p1 != NULL)
				for ( ; *p2 != 0 && *p2 != '~'; )
					p2++;

			bool bFailed = false;

			if ( wcslen( p1) != wcslen( p2) )
			{
				bFailed = true;
			} 
			else if ( p1 != NULL && p2 != NULL && *p1 != 0)
			{
				// skip starting ~
				p2++;
				for ( p1++ ; *p1 != 0; p1++)
				{
					const wchar_t * p = p2;

					for ( ; *p != 0 && *p != *p1; )
						p++;

					if ( *p == 0)
					{
						bFailed = true;
						break;
					}
				}
			}

			if ( bFailed )
			{
				wprintf(L"differences in menu item (~): %s (menu %s)\n", pItem1->GetID(), menu1->GetName());
				nRet = 1;
			}

			// find stringtable item for menu item
			string_map::iterator ii = strMap.find( pItem1->GetID() );
			if ( ii == strMap.end())
			{
				wprintf(L"missing string in stringtable for menu item: %s (menu %s)\n", pItem1->GetID(), menu1->GetName());
				nRet = 1;
			}
			else
			{
				std::wstring str = ii->second.substr(1, ii->second.length() -2); // trim quotes

				// string shall contain two blocks separated by new line '\n'
				// each block shall be more than one character long and not of the same characters
				// escape characters are not decoded
				p1 = str.c_str();

				if ( wcslen( p1) > 1)
				{
					bool b1 = true, bN = true, b2 = true;
					p2 = p1;
					for ( p1++; *p1 != 0; p1++)
					{
						if ( *p1 != *p2 )
						{
							if ( *p1 == '\\' && *++p1 == 'n')
							{
								p1++;
								bN = false;
								break;
							}
							b1 = false;
						}
					}

					if ( *p1 != 0)
						p2 = p1++;

					for ( ; *p1 != 0; p1++)
					{
						if ( *p1 != *p2)
						{
							b2 = false;
							break;
						}
					}

					bFailed = b1 || bN || b2;
				}
				else
					bFailed = true;

				if ( bFailed)
				{
					wprintf(L"improper string in stringtable for menu item: %s (menu %s)\n", pItem1->GetID(), menu1->GetName());
					nRet = 1;
				}
			}
		}
	}

	return nRet;
}

int compareMenus(RCFile_C *rc1, RCFile_C *rc2, const wchar_t * id)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	string_map str2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Menu_C), &pRes)) != -1)
		if ( ! wcscmp( id, pRes->GetName()) )
			map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Menu_C), &pRes)) != -1)
		if ( ! wcscmp( id, pRes->GetName()) )
			map2[pRes->GetName()] = pRes;

	// collect resources from rc2 - stringtable
	pRes = NULL;
	rc2->GetResource(0, &typeid(StringTable_C), &pRes);
	if (pRes != NULL) {
		StringTable_C *pTable = (StringTable_C*)pRes;
		String_C *pString;
		resIdx = 0;
		while ((pString = pTable->GetString(resIdx++)) != NULL)
			str2[pString->GetID()] = pString->GetText();
	}

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: MENU %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: MENU %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// compare menus
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) != map2.end()) {
			Menu_C *menu1 = (Menu_C*)it->second;
			Menu_C *menu2 = (Menu_C*)map2.find(it->first)->second;
			ret |= compareMenu(menu1, menu2, str2);
		}
	}

	return ret;
}

int compareFiles(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;

	ret |= compareMenus(rc1, rc2, L"IDR_HOTKEY");

	return ret;
}
