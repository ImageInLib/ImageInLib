#include "headers.h"
#include "merge.h"
#include "rcfile.h"
#include "typeinfo.h"
#include <mbstring.h>


extern bool g_verbose;
extern std::wstring g_mark;

wchar_t *szSourceRC = NULL;
size_t nSourceRCSize = 0;
size_t nLastTextReplacePos = 0;


std::wstring markText(const wchar_t *pszText, bool markEmpty = false) 
{
	if (pszText == NULL || *pszText == 0)
		return L"";
	if (wcscmp(pszText, L"\"\"") == 0 && !markEmpty)
		return L"\"\"";
	std::wstring text;
	if (*pszText == '"')
		text = L"\"" + g_mark + (pszText + 1);
	else
		text = pszText;
	return text;
}

void replaceText(FILE *pFile, size_t nReplaceStart, size_t nReplaceEnd, const wchar_t* pszNewText)
{
	if (nLastTextReplacePos < nReplaceStart)
		fwrite(szSourceRC + nLastTextReplacePos, 2, nReplaceStart - nLastTextReplacePos, pFile);
	if (*pszNewText)
		fwrite(pszNewText, 2, wcslen(pszNewText), pFile);
	nLastTextReplacePos = nReplaceEnd;
	fflush(pFile);
}

void mergeTextsFromDialog(Dialog_C *pDlg, RCFile_C *pTranslatedRC, FILE *pFile)
{
	Dialog_C *pNewDlg;
	int ret = 0;
	Control_C *pCtrl;
	int ctrlIdx;
	int newCtrlIdx;
	bool skip = false;

	pNewDlg = pTranslatedRC->FindDialog(pDlg->GetName());
	if (pNewDlg == NULL)
	{
		fwprintf(stdout, L"WARNING : dialog %s translation is missing\n", pDlg->GetName());
		skip = true;
	}

	if (skip)
	{
		replaceText(pFile, pDlg->GetCaptionStartPos(), pDlg->GetCaptionEndPos(), markText(pDlg->GetCaption()).c_str());
	}
	else
	{
		if (pNewDlg->GetTranslationHint() && _wcsicmp(pDlg->GetCaption(), pNewDlg->GetTranslationHint()) != 0)
		{
			fwprintf(stdout, L"WARNING : original caption of dialog %s is different. Dialog caption must be merged manually.\n", 
				pDlg->GetName());
			replaceText(pFile, pDlg->GetCaptionStartPos(), pDlg->GetCaptionEndPos(), markText(pDlg->GetCaption()).c_str());
		} else {
			replaceText(pFile, pDlg->GetCaptionStartPos(), pDlg->GetCaptionEndPos(), pNewDlg->GetCaption());
		}
	}

	// iterate controls
	ctrlIdx = 0;
	newCtrlIdx = 0;
	while ((pCtrl = pDlg->GetControl(ctrlIdx++)) != NULL)
	{
		if (skip) {
			replaceText(pFile, pCtrl->GetTextStartPos(), pCtrl->GetTextEndPos(), markText(pCtrl->GetText()).c_str());
			continue;
		}
		Control_C *pNewCtrl = pNewDlg->GetControl(newCtrlIdx);
		if (pNewCtrl == NULL) {
			fwprintf(stdout, L"WARNING : differences found in dialog %s. Dialog must be merged manually.\n", pDlg->GetName());
			replaceText(pFile, pCtrl->GetTextStartPos(), pCtrl->GetTextEndPos(), markText(pCtrl->GetText()).c_str());
			skip = true;
			continue;
		}
		if (wcslen(pCtrl->GetText()) == 0 || _wcsicmp(pCtrl->GetText(), L"\"\"") == 0)
		{
			newCtrlIdx++;
			continue;
		}
		if (_wcsicmp(pNewCtrl->GetType(), pCtrl->GetType()) != 0 || _wcsicmp(pNewCtrl->GetID(), pCtrl->GetID()) != 0) {
			fwprintf(stdout, L"WARNING : differences found in dialog %s. Dialog must be merged manually.\n", pDlg->GetName());
			replaceText(pFile, pCtrl->GetTextStartPos(), pCtrl->GetTextEndPos(), markText(pCtrl->GetText()).c_str());
			skip = true;
			continue;
		}
		if (pNewCtrl->GetTranslationHint() && _wcsicmp(pCtrl->GetText(), pNewCtrl->GetTranslationHint()) != 0) {
			fwprintf(stdout, L"WARNING : original text of control %s in dialog %s is different. Dialog control must be merged manually.\n", 
				pCtrl->GetID(), pDlg->GetName());
			replaceText(pFile, pCtrl->GetTextStartPos(), pCtrl->GetTextEndPos(), markText(pCtrl->GetText()).c_str());
		} else {
			replaceText(pFile, pCtrl->GetTextStartPos(), pCtrl->GetTextEndPos(), pNewCtrl->GetText());
		}
		newCtrlIdx++;
	}

	return;
}

void mergeTextsFromMenu(Menu_C *pMenu, RCFile_C *pTranslatedRC, FILE *pFile)
{
	Menu_C *pNewMenu;
	MenuPopup_C *pNewPopup = NULL;
	MenuPopup_C *pPopup;
	std::stack<MenuPopup_C*> newPopupStack;
	std::stack<int> popupStack;
	std::stack<int> popupIdxStack;
	int itemIdx;
	int popupID = 0;

	pNewMenu = pTranslatedRC->FindMenu(pMenu->GetName());

	// iterate items
	pNewPopup = NULL;
	pPopup = NULL;
	itemIdx = 0;
	for (;;)
	{
		MenuItem_C *pNewItem;
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
			popupID = popupIdxStack.top();
			popupIdxStack.pop();
			pPopup = pPopup->GetOwner();
			pNewPopup = newPopupStack.top();
			newPopupStack.pop();
			continue;
		}
		if (dynamic_cast<MenuPopup_C*>(pItem)) 
		{
			newPopupStack.push(pNewPopup);
			popupStack.push(itemIdx);
			itemIdx = 0;
			popupID++;
			popupIdxStack.push(popupID);
			wchar_t num[16];
			swprintf_s(num, _countof(num), L":#%d", popupID);
			pNewItem = NULL;
//			pNewPopup = NULL;
			pPopup = dynamic_cast<MenuPopup_C*>(pItem);
			if (pNewPopup == (MenuPopup_C*)1 || pNewMenu == NULL)
			{
				;
			}
			else if (pNewPopup)
			{
				for (int i = 0, newPopupID = 0; i < pNewPopup->GetItemsCount(); i++) {
					pNewItem = pNewPopup->GetItem(i);
					if (pNewItem == NULL)
						continue;
					if (dynamic_cast<MenuPopup_C*>(pNewItem) != NULL)
					{
						newPopupID++;
						if (pNewItem->GetTranslationHint())
						{
							if (_wcsicmp(pItem->GetText(), pNewItem->GetTranslationHint()) == 0)
								break;
						}
						else
						{
							if (newPopupID == popupID)
								break;
						}
					}
					pNewItem = NULL;
				}
				pNewPopup = dynamic_cast<MenuPopup_C*>(pNewItem);
				popupID = 0;
			}
			else
			{
				for (int i = 0, newPopupID = 0; i < pNewMenu->GetItemsCount(); i++)
				{
					pNewItem = pNewMenu->GetItem(i);
					if (pNewItem == NULL)
						continue;
					if (dynamic_cast<MenuPopup_C*>(pNewItem) != NULL)
					{
						newPopupID++;
						if (pNewItem->GetTranslationHint())
						{
							if (_wcsicmp(pItem->GetText(), pNewItem->GetTranslationHint()) == 0)
								break;
						}
						else
						{
							if (newPopupID == popupID)
								break;
						}
					}
					pNewItem = NULL;
				}
				pNewPopup = dynamic_cast<MenuPopup_C*>(pNewItem);
				popupID = 0;
			}
			if (pNewPopup == NULL)
			{
				pNewPopup = (MenuPopup_C*)1;
				replaceText(pFile, pPopup->GetTextStartPos(), pPopup->GetTextEndPos(), markText(pPopup->GetText()).c_str());
			}
			else if (pNewPopup != (MenuPopup_C*)1)
			{
				if (pNewPopup->GetTranslationHint() && _wcsicmp(pPopup->GetText(), pNewPopup->GetTranslationHint()) != 0)
				{
					fwprintf(stdout, L"WARNING : original text of menu popup %s in menu %s is different. Menu popup must be merged manually.\n", 
						pPopup->GetText(), pMenu->GetName());
					replaceText(pFile, pPopup->GetTextStartPos(), pPopup->GetTextEndPos(), markText(pPopup->GetText()).c_str());
				} else {
					replaceText(pFile, pPopup->GetTextStartPos(), pPopup->GetTextEndPos(), pNewPopup->GetText());
				}
			}
		} else if (_wcsicmp(pItem->GetID(), L"SEPARATOR") != 0) {
			pNewItem = NULL;
			if (pNewPopup == (MenuPopup_C*)1 || pNewMenu == NULL) {
				;
			} else if (pNewPopup) {
				for (int i = 0; i < pNewPopup->GetItemsCount(); i++) {
					pNewItem = pNewPopup->GetItem(i);
					if (pNewItem == NULL)
						continue;
					if (dynamic_cast<MenuPopup_C*>(pNewItem) == NULL && _wcsicmp(pItem->GetID(), pNewItem->GetID()) == 0)
						break;
					pNewItem = NULL;
				}
			} else {
				for (int i = 0; i < pNewMenu->GetItemsCount(); i++) {
					pNewItem = pNewMenu->GetItem(i);
					if (pNewItem == NULL)
						continue;
					if (dynamic_cast<MenuPopup_C*>(pNewItem) == NULL && _wcsicmp(pItem->GetID(), pNewItem->GetID()) == 0)
						break;
					pNewItem = NULL;
				}
			}
			if (pNewItem == NULL) {
				replaceText(pFile, pItem->GetTextStartPos(), pItem->GetTextEndPos(), markText(pItem->GetText()).c_str());
			} else {
				if (pNewItem->GetTranslationHint() && _wcsicmp(pItem->GetText(), pNewItem->GetTranslationHint()) != 0) {
					fwprintf(stdout, L"WARNING : original text of menu item %s in menu %s is different. Menu item must be merged manually.\n", 
						pItem->GetID(), pMenu->GetName());
					replaceText(pFile, pItem->GetTextStartPos(), pItem->GetTextEndPos(), markText(pItem->GetText()).c_str());
				} else {
					replaceText(pFile, pItem->GetTextStartPos(), pItem->GetTextEndPos(), pNewItem->GetText());
				}
			}
		}
	}

	return;
}

void mergeTextsFromStringTable(StringTable_C *pTable, RCFile_C *pTranslatedRC, FILE *pFile)
{
	int strIdx;

	String_C *pString;
	strIdx = 0;
	while ((pString = pTable->GetString(strIdx++)) != NULL) {
		String_C *pNewString = pTranslatedRC->FindString(pString->GetID());
		if (pNewString != NULL) {
			if (pNewString->GetTranslationHint() && _wcsicmp(pString->GetText(), pNewString->GetTranslationHint()) != 0) {
				fwprintf(stdout, L"WARNING : original text of string %s is different. String must be merged manually.\n", pString->GetID());
				replaceText(pFile, pString->GetTextStartPos(), pString->GetTextEndPos(), markText(pString->GetText()).c_str());
			} else {
				replaceText(pFile, pString->GetTextStartPos(), pString->GetTextEndPos(), pNewString->GetText());
			}
		} else {
			if (g_verbose) {
				fflush(stdout);
				fwprintf(stdout, L"NOTICE : string %s translation is missing\n", pString->GetID());
			}
			replaceText(pFile, pString->GetTextStartPos(), pString->GetTextEndPos(), markText(pString->GetText()).c_str());
		}
	}
}

int mergeTexts(RCFile_C *pSourceRC, RCFile_C *pTranslatedRC, const wchar_t *pszSourceRCFile, const wchar_t *pszOutputRCFile)
{
	Resource_C *pRes;
	int resIdx;

	// open source file
	FILE *in = NULL;
	if (_wfopen_s(&in, pszSourceRCFile, L"rb") != 0)
	{
		fwprintf(stderr, L"fatal error: Cannot open source file: '%s'\n", pszSourceRCFile);
		return NULL;
	}

	// detect file character set
	bool bUnicode = false;
	char c1, c2;
	c1 = fgetc(in);
	c2 = fgetc(in);
	fseek(in, 0, SEEK_SET);
	size_t len = _filelength(_fileno(in));
	if (len >= 2 && c1 == '\xff' && c2 == '\xfe') 
	{
		// Unicode with BOM
		bUnicode = true;
		len--;
	}
	else if (len >= 2 && c1 != 0 && c2 == 0) 
	{
		// Unicode without BOM
		bUnicode = true;
	}
	else
	{
		// 8-bit
		bUnicode = false;
	}

	// allocate parser buffer
	nSourceRCSize = (bUnicode ? len / sizeof(wchar_t) : len) + 1;
	szSourceRC = new wchar_t[nSourceRCSize];
	if (!szSourceRC)
	{
		fclose(in);
		fwprintf(stderr, L"fatal error: Cannot allocate buffer for source file: '%s'\n", pszSourceRCFile);
		return NULL;
	}
	memset(szSourceRC, 0, nSourceRCSize * sizeof(wchar_t));

	// read .RC file
	if (bUnicode)
	{
		for (size_t pos = 0, idx = 0; pos < len / sizeof(wchar_t); pos++)
		{
			wchar_t c = fgetwc(in);
			if (c == L'\xfeff')
			{
				nSourceRCSize--;
			}
			else if (c == '\r')
			{
				if (idx < len / sizeof(wchar_t) - 1)
				{
					c = fgetwc(in);
					if (c != '\n')
						szSourceRC[idx++] = '\r';
					else
						nSourceRCSize--;
					pos++;
				}
				szSourceRC[idx++] = c;
			}
			else
			{
				szSourceRC[idx++] = c;
			}
		}
	}
	else
	{
		for (size_t pos = 0, idx = 0; pos < len; pos++, idx++)
		{
			wchar_t c = _mbbtombc(fgetc(in));
			if (c == '\r')
			{
				if (idx < len / sizeof(wchar_t) - 1)
				{
					c = _mbbtombc(fgetc(in));
					if (c != '\n')
						szSourceRC[idx++] = '\r';
					else
						nSourceRCSize--;
					pos++;
				}
			}
			szSourceRC[idx] = c;
		}
	}

	// close source file
	fclose(in);
	in = NULL;

	// open output merged file
	FILE *pFile;
	if (_wfopen_s(&pFile, pszOutputRCFile, L"w, ccs=UTF-16LE") || pFile == NULL)
		return 1;
	if (pFile == NULL)
	{
		delete szSourceRC;
		return 1;
	}

	// iterate all resources
	for (resIdx = 0; resIdx < pSourceRC->GetResourceCount(); resIdx++)
	{
		pRes = pSourceRC->GetResource(resIdx);

		if (dynamic_cast<StringTable_C*>(pRes))
			mergeTextsFromStringTable(dynamic_cast<StringTable_C*>(pRes), pTranslatedRC, pFile);
		else if (dynamic_cast<Dialog_C*>(pRes))
			mergeTextsFromDialog(dynamic_cast<Dialog_C*>(pRes), pTranslatedRC, pFile);
		else if (dynamic_cast<Menu_C*>(pRes))
			mergeTextsFromMenu(dynamic_cast<Menu_C*>(pRes), pTranslatedRC, pFile);
	}
	replaceText(pFile, nSourceRCSize, nSourceRCSize, L"");

	delete szSourceRC;
	fclose(pFile);
	return 0;
}
