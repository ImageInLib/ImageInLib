#include "headers.h"
#include "compare.h"
#include "rcfile.h"
#include "typeinfo.h"

typedef std::map<std::wstring, Resource_C*> resources_map;
typedef std::map<std::wstring, std::wstring> string_map;
typedef std::list<std::wstring> string_list;

int compareCursors(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Cursor_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Cursor_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: CURSOR %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: CURSOR %s\n", it->second->GetName());
			ret = 1;
		}
	}

	return ret;
}

int compareBitmaps(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Bitmap_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Bitmap_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: BITMAP %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: BITMAP %s\n", it->second->GetName());
			ret = 1;
		}
	}

	return ret;
}

int compareControlOptions(const wchar_t *pOptions1, const wchar_t *pOptions2)
{
	if (wcscmp(pOptions1, pOptions2) == 0)
		return 0;

	while (*pOptions1 && *pOptions2)
	{
		while (*pOptions1 == ' ' || *pOptions1 == '\t')
			pOptions1++;
		while (*pOptions2 == ' ' || *pOptions2 == '\t')
			pOptions2++;

		// ignore numbers
		if (*pOptions1 == '-' || (*pOptions1 >= '0' && *pOptions1 <= '9'))
		{
			while (*pOptions1 == '-' || (*pOptions1 >= '0' && *pOptions1 <= '9'))
				pOptions1++;
			while (*pOptions2 == '-' || (*pOptions2 >= '0' && *pOptions2 <= '9'))
				pOptions2++;
			while (*pOptions1 == ' ' || *pOptions1 == '\t')
				pOptions1++;
			while (*pOptions2 == ' ' || *pOptions2 == '\t')
				pOptions2++;
		}

		if (*pOptions1 == ',')
		{
			while (*pOptions2 == ' ' || *pOptions2 == '\t')
				pOptions2++;
			if (*pOptions1 != *pOptions2)
			{
				return 1;
			}
			pOptions1++;
			pOptions2++;
		}
		else
		{
			while (*pOptions1 && *pOptions1 != ',')
			{
				while (*pOptions1 == ' ' || *pOptions1 == '\t')
					pOptions1++;
				while (*pOptions2 == ' ' || *pOptions2 == '\t')
					pOptions2++;
				if (*pOptions1 != *pOptions2)
				{
					return 1;
				}
				pOptions1++;
				pOptions2++;
			}
		}
	}

	while (*pOptions1 == ' ' || *pOptions1 == '\t')
		pOptions1++;
	while (*pOptions2 == ' ' || *pOptions2 == '\t')
		pOptions1++;

	if (*pOptions1 || *pOptions2)
		return 1;

	return 0;
}

int compareDialogStyle(const wchar_t *pOptions1, const wchar_t *pOptions2)
{
	if (wcscmp(pOptions1, pOptions2) == 0)
		return 0;

	std::set<std::wstring> styles1;
	std::set<std::wstring> styles2;

	while (*pOptions1)
	{
		// skip white spaces
		while (*pOptions1 == ' ' || *pOptions1 == '\t')
			pOptions1++;

		std::wstring style1;
		while (*pOptions1 && *pOptions1 != ' ' && *pOptions1 != '\t')
			style1 += *pOptions1++;
		if (!style1.empty())
			styles1.insert(style1);

		while (*pOptions1 == ' ' || *pOptions1 == '\t')
			pOptions1++;

		if (*pOptions1 != '|')
			break;
		pOptions1++;
	}
	while (*pOptions2)
	{
		// skip white spaces
		while (*pOptions2 == ' ' || *pOptions2 == '\t')
			pOptions2++;

		std::wstring style2;
		while (*pOptions2 && *pOptions2 != ' ' && *pOptions2 != '\t')
			style2 += *pOptions2++;
		if (!style2.empty())
			styles2.insert(style2);

		while (*pOptions2 == ' ' || *pOptions2 == '\t')
			pOptions2++;

		if (*pOptions2 != '|')
			break;
		pOptions2++;
	}

	std::set<std::wstring>::iterator it;
	for (it = styles1.begin(); it != styles1.end(); it++)
	{
		if (styles2.find(*it) == styles2.end())
			return 1;
	}
	for (it = styles2.begin(); it != styles2.end(); it++)
	{
		if (styles1.find(*it) == styles1.end())
			return 1;
	}

	return 0;
}

std::wstring getDialogStyle(const wchar_t *pOptions)
{
	std::wstring style;

	pOptions = wcsstr(pOptions, L"STYLE");
	if (pOptions != NULL)
		pOptions += wcslen(L"STYLE");
	if (!pOptions)
		return style;

	while (*pOptions)
	{
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		while (*pOptions && *pOptions != ' ' && *pOptions != '\t')
			style += *pOptions++;
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		if (*pOptions != '|')
			break;
		pOptions++;
		style += L" | ";
	}

	return style;
}

std::wstring getDialogExStyle(const wchar_t *pOptions)
{
	std::wstring style;

	pOptions = wcsstr(pOptions, L" EXSTYLE ");
	if (pOptions != NULL)
		pOptions += wcslen(L" EXSTYLE ");
	if (!pOptions)
		return style;

	while (*pOptions)
	{
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		while (*pOptions && *pOptions != ' ' && *pOptions != '\t')
			style += *pOptions++;
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		if (*pOptions != '|')
			break;
		pOptions++;
		style += L" | ";
	}

	return style;
}

std::wstring getDialogFont(const wchar_t *pOptions)
{
	std::wstring font;

	pOptions = wcsstr(pOptions, L" FONT ");
	if (pOptions != NULL)
		pOptions += wcslen(L" FONT ");
	if (!pOptions)
		return font;

	while (*pOptions)
	{
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		if (*pOptions == '"')
		{
			font += *pOptions++;
			while (*pOptions && *pOptions != '"')
				font += *pOptions++;
			if (*pOptions)
				font += *pOptions++;
		}
		else
		{
			while (*pOptions && *pOptions != ' ' && *pOptions != '\t')
				font += *pOptions++;
		}
		while (*pOptions == ' ' || *pOptions == '\t')
			pOptions++;
		if (*pOptions != ',')
			break;
		pOptions++;
		font += L", ";
	}

	return font;
}

int compareDialog(Dialog_C *dlg1, Dialog_C *dlg2)
{
	int ret = 0;
	Control_C *pCtrl1, *pCtrl2;
	typedef std::map<std::wstring, Control_C*> controls_map;
	controls_map::iterator it;
	controls_map map1, map2;
	int ctrlIdx;
	int staticIdx;
	wchar_t buf[32];

	if (wcscmp(dlg1->GetType(), dlg2->GetType()) != 0)
	{
		wprintf(L"different dialog types in %s ('%s' vs. '%s')\n", dlg1->GetName(), dlg1->GetType(), dlg2->GetType());
		ret = 1;
	}

	std::wstring dlgStyle1 = getDialogStyle(dlg1->GetOptions());
	std::wstring dlgStyle2 = getDialogStyle(dlg2->GetOptions());
	if (compareDialogStyle(dlgStyle1.c_str(), dlgStyle2.c_str()) != 0)
	{
		wprintf(L"different dialog styles in %s ('%s' vs. '%s')\n", dlg1->GetName(), dlgStyle1.c_str(), dlgStyle2.c_str());
		ret = 1;
	}

	std::wstring dlgExStyle1 = getDialogExStyle(dlg1->GetOptions());
	std::wstring dlgExStyle2 = getDialogExStyle(dlg2->GetOptions());
	if (compareDialogStyle(dlgExStyle1.c_str(), dlgExStyle2.c_str()) != 0)
	{
		wprintf(L"different dialog extended styles in %s ('%s' vs. '%s')\n", dlg1->GetName(), dlgExStyle1.c_str(), dlgExStyle2.c_str());
		ret = 1;
	}

	std::wstring dlgFont1 = getDialogFont(dlg1->GetOptions());
	std::wstring dlgFont2 = getDialogFont(dlg2->GetOptions());
	if (wcscmp(dlgFont1.c_str(), dlgFont2.c_str()) != 0)
	{
		wprintf(L"different dialog fonts in %s ('%s' vs. '%s')\n", dlg1->GetName(), dlgFont1.c_str(), dlgFont2.c_str());
		ret = 1;
	}

	if (ret != 0)
		return ret;

	// collect controls from dlg1
	ctrlIdx = 0;
	staticIdx = 0;
	while ((pCtrl1 = dlg1->GetControl(ctrlIdx++)) != NULL) 
	{
		controls_map::iterator it = map1.find(pCtrl1->GetID());
		if (it != map1.end())
			swprintf_s(buf, _countof(buf), L"#%d", staticIdx++);
		else
			buf[0] = '\0';
		map1[std::wstring(pCtrl1->GetID()) + buf] = pCtrl1;
	}

	// collect controls from dlg2
	ctrlIdx = 0;
	staticIdx = 0;
	while ((pCtrl2 = dlg2->GetControl(ctrlIdx++)) != NULL) 
	{
		controls_map::iterator it = map2.find(pCtrl2->GetID());
		if (it != map2.end())
			swprintf_s(buf, _countof(buf), L"#%d", staticIdx++);
		else
			buf[0] = '\0';
		map2[std::wstring(pCtrl2->GetID()) + buf] = pCtrl2;
	}

	// find missing controls
	for (it = map1.begin(); it != map1.end(); it++) 
	{
		if (map2.find(it->first) == map2.end()) 
		{
			wprintf(L"missing dialog control: %s %s in %s\n", it->second->GetType(), it->second->GetID(), dlg1->GetName());
			ret = 1;
		}
	}

	// find extra controls
	for (it = map2.begin(); it != map2.end(); it++) 
	{
		if (map1.find(it->first) == map1.end()) 
		{
			wprintf(L"extra dialog control: %s %s in %s\n", it->second->GetType(), it->second->GetID(), dlg1->GetName());
			ret = 1;
		}
	}

	if (ret == 0) {
		if (dlg1->GetControlsCount() < dlg2->GetControlsCount())
		{
			wprintf(L"missing dialog controls: %d controls in %s\n", dlg2->GetControlsCount() - dlg1->GetControlsCount(), dlg1->GetName());
			ret = 1;
		} 
		else if (dlg1->GetControlsCount() > dlg2->GetControlsCount()) 
		{
			wprintf(L"extra dialog controls: %d controls in %s\n", dlg1->GetControlsCount() - dlg2->GetControlsCount(), dlg1->GetName());
			ret = 1;
		}
	}

	// compare control types
	for (it = map1.begin(); it != map1.end(); it++) 
	{
		if (map2.find(it->first) != map2.end()) 
		{
			Control_C *pCtrl1 = (Control_C*)it->second;
			Control_C *pCtrl2 = (Control_C*)map2.find(it->first)->second;
			if (wcscmp(pCtrl1->GetType(), pCtrl2->GetType()) != 0) 
			{
				wprintf(L"different control types ('%s' vs. '%s'): %s in %s\n", pCtrl1->GetType(), pCtrl2->GetType(), pCtrl1->GetID(), dlg1->GetName());
				ret = 1;
			}
			if (compareControlOptions(pCtrl1->GetOptions(), pCtrl2->GetOptions()) != 0)
			{
				wprintf(L"different control options ('%s' vs. '%s'): %s in %s\n", pCtrl1->GetOptions(), pCtrl2->GetOptions(), pCtrl1->GetID(), dlg1->GetName());
				ret = 1;
			}
		}
	}

	// compare control types
	if (ret == 0) 
	{
		for (ctrlIdx = 0; 
			 NULL != (pCtrl1 = dlg1->GetControl(ctrlIdx))&&
			 NULL != (pCtrl2 = dlg2->GetControl(ctrlIdx));
			 ctrlIdx++) 
		{
			if (wcscmp(pCtrl1->GetType(), pCtrl2->GetType()) != 0 ||
				wcscmp(pCtrl1->GetID(), pCtrl2->GetID()) != 0)
			{
				wprintf(L"different tab order in dialog %s\n", dlg1->GetName());
				ret = 1;
				break;
			}
		}
	}

	return ret;
}

int compareDialogs(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Dialog_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Dialog_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: DIALOG %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: DIALOG %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// compare dialogs
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) != map2.end()) {
			Dialog_C *dlg1 = (Dialog_C*)it->second;
			Dialog_C *dlg2 = (Dialog_C*)map2.find(it->first)->second;
			ret |= compareDialog(dlg1, dlg2);
		}
	}

	return ret;
}

int compareFonts(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Font_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Font_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: FONT %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: FONT %s\n", it->second->GetName());
			ret = 1;
		}
	}

	return ret;
}

int compareIcons(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Icon_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Icon_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: ICON %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: ICON %s\n", it->second->GetName());
			ret = 1;
		}
	}

	return ret;
}

int compareMenu(Menu_C *menu1, Menu_C *menu2)
{
	MenuPopup_C *pPopup1, *pPopup2;
	stack<int> popupStack;
	int itemIdx;

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
	}

	return 0;
}

int compareMenus(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Menu_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Menu_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

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
			ret |= compareMenu(menu1, menu2);
		}
	}

	return ret;
}

int compareRCDatas(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(RCData_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(RCData_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: RCDATA %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: RCDATA %s\n", it->second->GetName());
			ret = 1;
		}
	}

	return ret;
}

int parseFormatSpecification(const wchar_t *pString, std::list<std::wstring> &formatDirectives)
{
	// parse format specification
	std::wstring formatDirective;
	const wchar_t* formatSpec = pString;

	int argumentsCnt = 0;
	int state = 0;
	for (;;) 
	{
		if (state == 0)
		{

			// ordinary string
			if (*formatSpec == 0)
			{
				break;
			}
			if (*formatSpec == '%')
			{
				formatDirective = L"%";
				state = 2;
			} 
			else if (*formatSpec == '\\')
			{
				formatDirective = L"";
				state = 1;
			} 
			else 
			{
				formatDirective += *formatSpec;
			}
			formatSpec++;

		} 
		else if (state == 1) 
		{

			// escape sequences
			if (*formatSpec == 0) 
			{
				if (formatDirective.size() == 0) 
				{
					break;
				}
				int chr;
				if (formatDirective[0] == 'x')
					swscanf_s(formatDirective.c_str(), L"x%x", &chr);
				else
					swscanf_s(formatDirective.c_str(), L"%o", &chr);
				wchar_t str[2];
				str[0] = (char)chr;
				str[1] = 0;
				break;
			}
			wchar_t str[2] = {0, 0};
			if (formatDirective.size() == 0)
			{
				switch(*formatSpec) 
				{
				case 'a': *str = '\a'; break;
				case 'b': *str = '\a'; break;
				case 'f': *str = '\a'; break;
				case 'n': *str = '\a'; break;
				case 'r': *str = '\a'; break;
				case 't': *str = '\a'; break;
				case 'v': *str = '\a'; break;
				case '0': case '1':	case '2':
				case '3': case '4':	case '5':
				case '6': case '7':	case 'x': 
					formatDirective += *formatSpec;	break;
				default:
					*str = *formatSpec; break;
				}
				formatSpec++;
			} 
			else if (formatDirective[0] == 'x') 
			{
				if ((*formatSpec >= '0' && *formatSpec <= '9') 
				 || (*formatSpec >= 'a' && *formatSpec <= 'f') 
				 || (*formatSpec >= 'F' && *formatSpec <= 'F'))
				{
					formatDirective += *formatSpec; 
					formatSpec++;
				} 
				else
				{
					int chr;
					swscanf_s(formatDirective.c_str(), L"x%x", &chr);
					*str = (char)chr;
				}
			}
			else
			{
				if (*formatSpec >= '0' && *formatSpec <= '7')
				{
					formatDirective += *formatSpec; 
					formatSpec++;
				}
				else
				{
					int chr;
					swscanf_s(formatDirective.c_str(), L"%o", &chr);
					*str = (char)chr;
				}
			}
			if (*str)
			{
				formatDirective = L"";
				state = 0;
			}

		}
		else if (state == 2)
		{

			// format flags
			if (*formatSpec == 0)
			{
				break;
			}
			if (*formatSpec == '%')
			{
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}
			else if (*formatSpec == '+' 
				    || *formatSpec == '-' 
				    || *formatSpec == ' ' 
					|| *formatSpec == '0'
					|| *formatSpec == '#') 
			{
				formatDirective += *formatSpec;
				formatSpec++;
			}
			else if (*formatSpec >= '1' && *formatSpec <= '9')
			{
				state = 3;
			}
			else if (*formatSpec == '.')
			{
				formatDirective += *formatSpec;
				formatSpec++;
				state = 4;
			}
			else if (*formatSpec == 'l')
			{
				formatSpec++;
				if (*formatSpec == 0)
				{
					break;
				}
				else if (*formatSpec == 'l')
				{
					// ll size
					formatSpec++;
					if (*formatSpec == 0)
					{
						break;
					}
					else if (*formatSpec == 'd'
							|| *formatSpec == 'i'
							|| *formatSpec == 'o'
							|| *formatSpec == 'u'
							|| *formatSpec == 'x'
							|| *formatSpec == 'X') 
					{
						formatDirective += L"ll";
						state = 5;
					}
					else
					{
						formatDirective = L"";
						formatSpec++;
						state = 0;
					}
				}
				else if (*formatSpec == 'c'
					|| *formatSpec == 'C' 
				    || *formatSpec == 'd'
				    || *formatSpec == 'i'
				    || *formatSpec == 'o'
				    || *formatSpec == 'u'
				    || *formatSpec == 'x'
				    || *formatSpec == 'X'
				    || *formatSpec == 'f' 
				    || *formatSpec == 's'
				    || *formatSpec == 'S') 
				{
					formatDirective += L"l";
					state = 5;
				}
				else
				{
					formatDirective = L"";
					formatSpec++;
					state = 0;
				}
			}
			else if (*formatSpec == 'I')
			{
				formatSpec++;
				if (*formatSpec == 0)
				{
					break;
				}
				else if ((*formatSpec == '3' && *(formatSpec+1) == '2')
					  || (*formatSpec == '6' && *(formatSpec+1) == '4'))
				{
					// I32 or I64 size
					formatDirective += *formatSpec++;
					formatDirective += *formatSpec++;
					if (*formatSpec == 0)
					{
						break;
					}
					else if (*formatSpec == 'd'
							|| *formatSpec == 'i'
							|| *formatSpec == 'o'
							|| *formatSpec == 'u'
							|| *formatSpec == 'x'
							|| *formatSpec == 'X') 
					{
						state = 5;
					}
					else
					{
						formatDirective = L"";
						formatSpec++;
						state = 0;
					}
				}
				else if (*formatSpec == 'd'
						|| *formatSpec == 'i'
						|| *formatSpec == 'o'
						|| *formatSpec == 'u'
						|| *formatSpec == 'x'
						|| *formatSpec == 'X') 
				{
					formatDirective += L"I";
					state = 5;
				}
				else
				{
					formatDirective = L"";
					formatSpec++;
					state = 0;
				}
			}
			else if (*formatSpec == 'h')
			{
				formatSpec++;
				if (*formatSpec == 0)
				{
					break;
				}
				else if (*formatSpec == 'c'
						|| *formatSpec == 'C' 
						|| *formatSpec == 'd'
						|| *formatSpec == 'i'
						|| *formatSpec == 'o'
						|| *formatSpec == 'u'
						|| *formatSpec == 'x'
						|| *formatSpec == 'X'
						|| *formatSpec == 's'
						|| *formatSpec == 'S') 
				{
					formatDirective += L"h";
					state = 5;
				}
				else
				{
					formatDirective = L"";
					formatSpec++;
					state = 0;
				}
			}
			else if (*formatSpec == 'L')
			{
				formatSpec++;
				if (*formatSpec == 0)
				{
					break;
				}
				else if (*formatSpec == 'f') 
				{
					formatDirective += L"L";
					state = 5;
				}
				else
				{
					formatDirective = L"";
					formatSpec++;
					state = 0;
				}
			}
			else if (*formatSpec == 'w')
			{
				formatSpec++;
				if (*formatSpec == 0)
				{
					break;
				}
				else if (*formatSpec == 'c' || *formatSpec == 's') 
				{
					formatDirective += L"w";
					state = 5;
				}
				else
				{
					formatDirective = L"";
					formatSpec++;
					state = 0;
				}
			}
			else if (*formatSpec == 'c' 
					|| *formatSpec == 'C' 
				    || *formatSpec == 'd'
				    || *formatSpec == 'i'
				    || *formatSpec == 'o'
				    || *formatSpec == 'u'
				    || *formatSpec == 'x'
				    || *formatSpec == 'X'
				    || *formatSpec == 'e'
				    || *formatSpec == 'E'
				    || *formatSpec == 'f'
				    || *formatSpec == 'F'
				    || *formatSpec == 'g'
				    || *formatSpec == 'G'
				    || *formatSpec == 's'
				    || *formatSpec == 'S') 
			{
				state = 5;
			}
			else
			{
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}

		}
		else if (state == 3)
		{

			// width specification
			if (*formatSpec == 0)
			{
				break;
			}
			if (*formatSpec >= '0' && *formatSpec <= '9')
			{
				formatDirective += *formatSpec;
				formatSpec++;
			}
			else if (*formatSpec == '.')
			{
				formatDirective += *formatSpec;
				formatSpec++;
				state = 4;
			}
			else if (*formatSpec == 'c' 
					|| *formatSpec == 'C' 
				    || *formatSpec == 'd'
				    || *formatSpec == 'i'
				    || *formatSpec == 'o'
				    || *formatSpec == 'u'
				    || *formatSpec == 'x'
				    || *formatSpec == 'X'
				    || *formatSpec == 'e'
				    || *formatSpec == 'E'
				    || *formatSpec == 'f'
				    || *formatSpec == 'F'
				    || *formatSpec == 'g'
				    || *formatSpec == 'G'
				    || *formatSpec == 's'
				    || *formatSpec == 'S') 
			{
				state = 5;
			}
			else
			{
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}

		}
		else if (state == 4)
		{

			// precision specification
			if (*formatSpec == 0)
			{
				break;
			}
			if (*formatSpec >= '0' && *formatSpec <= '9')
			{
				formatDirective += *formatSpec;
				formatSpec++;
			}
			else if (*formatSpec == 'c' 
					|| *formatSpec == 'C' 
				    || *formatSpec == 'd'
				    || *formatSpec == 'i'
				    || *formatSpec == 'o'
				    || *formatSpec == 'u'
				    || *formatSpec == 'x'
				    || *formatSpec == 'X'
				    || *formatSpec == 'e'
				    || *formatSpec == 'E'
				    || *formatSpec == 'f'
				    || *formatSpec == 'F'
				    || *formatSpec == 'g'
				    || *formatSpec == 'G'
				    || *formatSpec == 's'
				    || *formatSpec == 'S') 
			{
				state = 5;
			}
			else
			{
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}

		}
		else if (state == 5)
		{

			// type specification
			if (*formatSpec == 0)
			{
				break;
			}
			if (*formatSpec == 'c' 
				|| *formatSpec == 'C' 
			    || *formatSpec == 'd'
			    || *formatSpec == 'i'
			    || *formatSpec == 'o'
			    || *formatSpec == 'u'
			    || *formatSpec == 'x'
			    || *formatSpec == 'X'
			    || *formatSpec == 'e'
			    || *formatSpec == 'E'
			    || *formatSpec == 'f'
			    || *formatSpec == 'F'
			    || *formatSpec == 'g'
			    || *formatSpec == 'G'
				|| *formatSpec == 's'
				|| *formatSpec == 'S') 
			{
				formatDirective += *formatSpec;
				formatDirectives.push_back(formatDirective);
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}
			else
			{
				formatDirective = L"";
				formatSpec++;
				state = 0;
			}
		}
	}

	return 0;
}

int compareStringTable(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	string_map::iterator it;
	string_map map1, map2;
	Resource_C *pRes;
	int strIdx;

	// collect resources from rc1
	pRes = NULL;
	rc1->GetResource(0, &typeid(StringTable_C), &pRes);
	if (pRes != NULL) {
		StringTable_C *pTable = (StringTable_C*)pRes;
		String_C *pString;
		strIdx = 0;
		while ((pString = pTable->GetString(strIdx++)) != NULL)
			map1[pString->GetID()] = pString->GetText();
	}

	// collect resources from rc2
	pRes = NULL;
	rc2->GetResource(0, &typeid(StringTable_C), &pRes);
	if (pRes != NULL) {
		StringTable_C *pTable = (StringTable_C*)pRes;
		String_C *pString;
		strIdx = 0;
		while ((pString = pTable->GetString(strIdx++)) != NULL)
			map2[pString->GetID()] = pString->GetText();
	}

	// find missing strings and compare string format specifiers
	for (it = map1.begin(); it != map1.end(); it++) {
		string_map::iterator it2 = map2.find(it->first);
		if (it2 == map2.end()) {
			wprintf(L"missing string: %s\n", it->first.c_str());
			ret = 1;
		}
		else
		{
			string_list formatDirectives1;
			string_list formatDirectives2;
			parseFormatSpecification(it->second.c_str(), formatDirectives1);
			parseFormatSpecification(it2->second.c_str(), formatDirectives2);
			if (formatDirectives1.size() != formatDirectives2.size())
			{
				wprintf(L"different string format: %s ('%s' vs. '%s')\n", 
					it->first.c_str(),
					it->second.c_str(),
					it2->second.c_str()
					);
				ret = 1;
			}
			else
			{
				for (string_list::iterator fit1 = formatDirectives1.begin(), fit2 = formatDirectives2.begin(); 
					fit1 != formatDirectives1.end() && fit2 != formatDirectives2.end(); 
					fit1++, fit2++)
				{
					std::wstring fmtSpec1 = *fit1;
					std::wstring fmtSpec2 = *fit2;
					if (*fit1 != *fit2)
					{
						wprintf(L"different string format: %s ('%s' vs. '%s')\n", 
							it->first.c_str(),
							it->second.c_str(),
							it2->second.c_str()
							);
						ret = 1;
						break;
					}
				}
			}
		}
	}

	// find extra strings
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra string: %s\n", it->first.c_str());
			ret = 1;
		}
	}

	return ret;
}

int compareToolbar(Toolbar_C *bar1, Toolbar_C *bar2)
{
	int itemIdx;

	if (bar1->GetItemsCount() != bar2->GetItemsCount()) {
		wprintf(L"differences in toolbar: %s\n", bar1->GetName());
		return 1;
	}

	// compare toolbar
	itemIdx = 0;
	for (;;) {
		const wchar_t *item1, *item2;
		item1 = bar1->GetItem(itemIdx);
		item2 = bar2->GetItem(itemIdx);
		itemIdx++;
		if (item1 == NULL)
			break;
		if (wcscmp(item1, item2) != 0) {
			wprintf(L"differences in toolbar: %s\n", bar1->GetName());
			return 1;
		}
	}

	return 0;
}

int compareToolbars(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(Toolbar_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(Toolbar_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: TOOLBAR %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: TOOLBAR %s\n", it->second->GetName());
			ret = 1;
		}
	}

	// compare toolbars
	for (it = map1.begin(); it != map1.end(); it++) {
		if (map2.find(it->first) != map2.end()) {
			Toolbar_C *bar1 = (Toolbar_C*)it->second;
			Toolbar_C *bar2 = (Toolbar_C*)map2.find(it->first)->second;
			ret |= compareToolbar(bar1, bar2);
		}
	}

	return ret;
}

int compareUserResources(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;
	resources_map::iterator it;
	resources_map map1, map2;
	Resource_C *pRes;
	int resIdx;

	// collect resources from rc1
	resIdx = 0;
	while ((resIdx = rc1->GetResource(resIdx, &typeid(UserDef_C), &pRes)) != -1)
		map1[pRes->GetName()] = pRes;

	// collect resources from rc2
	resIdx = 0;
	while ((resIdx = rc2->GetResource(resIdx, &typeid(UserDef_C), &pRes)) != -1)
		map2[pRes->GetName()] = pRes;

	// find missing resources
	for (it = map1.begin(); it != map1.end(); it++) {
		UserDef_C *pRes = (UserDef_C*)it->second;
		if (map2.find(it->first) == map2.end()) {
			wprintf(L"missing resource: %s %s\n", pRes->GetType(), pRes->GetName());
			ret = 1;
		}
	}

	// find extra resources
	for (it = map2.begin(); it != map2.end(); it++) {
		UserDef_C *pRes = (UserDef_C*)it->second;
		if (map1.find(it->first) == map1.end()) {
			wprintf(L"extra resource: %s %s\n", pRes->GetType(), pRes->GetName());
			ret = 1;
		}
	}

	return ret;
}

int compareFiles(RCFile_C *rc1, RCFile_C *rc2)
{
	int ret = 0;

	ret |= compareBitmaps(rc1, rc2);
	ret |= compareCursors(rc1, rc2);
	ret |= compareDialogs(rc1, rc2);
	ret |= compareFonts(rc1, rc2);
	ret |= compareIcons(rc1, rc2);
	ret |= compareMenus(rc1, rc2);
	ret |= compareRCDatas(rc1, rc2);
	ret |= compareStringTable(rc1, rc2);
	ret |= compareToolbars(rc1, rc2);
	ret |= compareUserResources(rc1, rc2);

	return ret;
}
