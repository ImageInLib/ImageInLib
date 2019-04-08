#include "headers.h"
#include "extract.h"
#include "rcfile.h"
#include "typeinfo.h"


std::wstring unquoteString(const wchar_t *str)
{
	std::wstring s;
	if (*str == '"')
	{
		str++;
		while (*str)
		{
			if (*str == '\\' && *(str+1) == '"') {
				s += L"\"";
				str++;
			} else if (*str == '"' && *(str+1) == '-') {
				s += L"'-";
				str++;
//			} else if (*str == '"' && *(str+1) >= '0' && *(str+1) <= '9' ) {
//				s += '\''; s+= *(str+1);
//				str++;
			} else if (*str == '"' && *(str+1) == '"') {
				s += L"\"";
				str++;
			} else if (*str == '"' && *(str+1) == 0) {
				;
			} else {
				s += *str;
			}
			str++;
		}
	}
	else
	{
		s = str;
	}
	return s;
}

int extractTextsFromDialog(Dialog_C *pDlg, FILE *pFile)
{
	int ret = 0;
	Control_C *pCtrl;
	std::wstring caption;
	int ctrlIdx;

	const wchar_t *p;
	p = wcsstr(pDlg->GetOptions(), L"CAPTION \"");
	if (p)
	{
		p += wcslen(L"CAPTION \"");
		while (*p != '"' && *p)
			caption += *p++;
	}

	fwprintf(pFile, L"DIALOG\t%s\t%s\n", pDlg->GetName(), caption.c_str());

	// iterate controls
	ctrlIdx = 0;
	while ((pCtrl = pDlg->GetControl(ctrlIdx++)) != NULL)
	{
		std::wstring text = unquoteString(pCtrl->GetText());
//		if (text != "" && text != "[HELMUT]")
			fwprintf(pFile, L"%s\t%s:%s\t%s\n", pCtrl->GetType(), pDlg->GetName(), pCtrl->GetID(), text.c_str());
	}

	return ret;
}

int extractTextsFromMenu(Menu_C *pMenu, FILE *pFile)
{
	MenuPopup_C *pPopup;
	std::stack<int> popupStack;
	std::stack<std::wstring> idStack;
	int itemIdx;
	int popupID = 0;

	idStack.push(pMenu->GetName());
	fwprintf(pFile, L"MENU\t%s\t\n", idStack.top().c_str());

	// iterate items
	pPopup = NULL;
	itemIdx = 0;
	for (;;)
	{
		MenuItem_C *pItem;
		if (pPopup)
			pItem = pPopup->GetItem(itemIdx);
		else
			pItem = pMenu->GetItem(itemIdx);
		itemIdx++;
		if (pItem == NULL)
		{
			if (popupStack.size() == 0)
				break;
			itemIdx = popupStack.top();
			popupStack.pop();
			idStack.pop();
			pPopup = pPopup->GetOwner();
			continue;
		}
		if (dynamic_cast<MenuPopup_C*>(pItem))
		{
			popupStack.push(itemIdx);
			pPopup = dynamic_cast<MenuPopup_C*>(pItem);
			itemIdx = 0;
			popupID++;
			wchar_t num[16];
			swprintf(num, _countof(num), L":#%d", popupID);
			std::wstring id = idStack.top() + num;
			idStack.push(id);
			fwprintf(pFile, L"POPUP\t%s\t%s\n", idStack.top().c_str(), unquoteString(pPopup->GetText()).c_str());
		}
		else
		{
			fwprintf(pFile, L"MENUITEM\t%s:%s\t%s\n", idStack.top().c_str(), pItem->GetID(), unquoteString(pItem->GetText()).c_str());
		}
	}

	return 0;
}

int extractTextsFromString(String_C *pString, FILE *pFile)
{
	fwprintf(pFile, L"STRING\t%s\t%s\n", pString->GetID(), unquoteString(pString->GetText()).c_str());

	return 0;
}

int extractTexts(RCFile_C *pRC, const wchar_t *pszTextFile)
{
	Resource_C *pRes;
	FILE *pFile;
	int resIdx;

	// open output text file
	if (_wfopen_s(&pFile, pszTextFile, L"w, ccs=UTF-8") || pFile == NULL)
		return 1;

	std::multimap<std::wstring, Menu_C*> menus;
	std::multimap<std::wstring, Dialog_C*> dialogs;
	std::multimap<std::wstring, String_C*> strings;

	for (resIdx = 0; resIdx < pRC->GetResourceCount(); resIdx++)
	{
		pRes = pRC->GetResource(resIdx);
		if (dynamic_cast<StringTable_C*>(pRes))
		{
			StringTable_C *pTable = dynamic_cast<StringTable_C*>(pRes);
			int strIdx = 0;
			String_C *pString;
			while ((pString = pTable->GetString(strIdx++)) != NULL)
			{
				std::pair<std::wstring, String_C*> p(pString->GetID(), pString);
				strings.insert(p);
			}
		}
		else if (dynamic_cast<Dialog_C*>(pRes))
		{
			Dialog_C* pDlg = dynamic_cast<Dialog_C*>(pRes);
			pair<std::wstring, Dialog_C*> p(pDlg->GetName(), pDlg);
			dialogs.insert(p);
		}
		else if (dynamic_cast<Menu_C*>(pRes))
		{
			Menu_C* pMenu = dynamic_cast<Menu_C*>(pRes);
			std::pair<std::wstring, Menu_C*> p(pMenu->GetName(), pMenu);
			menus.insert(p);
		}
	}

	std::multimap<std::wstring, Menu_C*>::iterator mit;
	for (mit = menus.begin(); mit != menus.end(); mit++)
		extractTextsFromMenu(mit->second, pFile);

	std::multimap<std::wstring, Dialog_C*>::iterator dit;
	for (dit = dialogs.begin(); dit != dialogs.end(); dit++)
		extractTextsFromDialog(dit->second, pFile);

	std::multimap<std::wstring, String_C*>::iterator sit;
	for (sit = strings.begin(); sit != strings.end(); sit++)
		extractTextsFromString(sit->second, pFile);

	return 0;
}
