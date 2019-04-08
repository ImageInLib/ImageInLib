// rcfile.cpp: implementation of the RCFile_C class.
//
//////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "rcfile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RCFile_C::RCFile_C()
{

}

RCFile_C::~RCFile_C()
{
	for (int i = 0; i < (int)m_resources.size(); i++)
		delete m_resources[i];
	m_resources.clear();
}

void RCFile_C::AddResource(Resource_C *pResource)
{
	m_resources.push_back(pResource);
	if (dynamic_cast<Menu_C*>(pResource) != NULL) {
		Menu_C *pMenu = dynamic_cast<Menu_C*>(pResource);
		m_menus[pMenu->GetName()] = pMenu;
	}
	if (dynamic_cast<Dialog_C*>(pResource) != NULL) {
		Dialog_C *pDialog = dynamic_cast<Dialog_C*>(pResource);
		m_dialogs[pDialog->GetName()] = pDialog;
	}
}

int RCFile_C::GetResourceCount()
{
	return m_resources.size();
}

Resource_C *RCFile_C::GetResource(int idx)
{
	if (idx < 0 || idx >= (int)m_resources.size())
		return NULL;
	else
		return m_resources[idx];
}

int RCFile_C::GetResource(int idx, const type_info *type, Resource_C **ppResource)
{
	if (idx < 0 || idx >= (int)m_resources.size()) {
		*ppResource = NULL;
		return -1;
	}
	while (idx < (int)m_resources.size()) {
		Resource_C *pResource;
		pResource = m_resources[idx++];
		if (typeid(*pResource) == *type) {
			*ppResource = pResource;
			return idx;
		}
	}
	if (idx >= (int)m_resources.size())
		idx = -1;
	return idx;
}


void RCFile_C::AddString(const wchar_t* stringID, String_C *pString)
{
	m_strings[stringID] = pString;
}


String_C* RCFile_C::FindString(const wchar_t *stringID)
{
	std::map<std::wstring, String_C*>::iterator it;
	it = m_strings.find(stringID);
	if (it != m_strings.end())
		return it->second;
	else
		return NULL;
}

Menu_C* RCFile_C::FindMenu(const wchar_t *menuID)
{
	std::map<std::wstring, Menu_C*>::iterator it;
	it = m_menus.find(menuID);
	if (it != m_menus.end())
		return it->second;
	else
		return NULL;
}

Dialog_C* RCFile_C::FindDialog(const wchar_t *dialogID)
{
	std::map<std::wstring, Dialog_C*>::iterator it;
	it = m_dialogs.find(dialogID);
	if (it != m_dialogs.end())
		return it->second;
	else
		return NULL;
}



//----------------------------------------
// Resource_C
//----------------------------------------
Resource_C::Resource_C(RCFile_C *pRCFile, const wchar_t *pszName) 
{ 
	m_pRCFile = pRCFile;
	m_sName = pszName; 
}

Resource_C::~Resource_C()
{

}

const wchar_t* Resource_C::GetName() 
{ 
	return m_sName.c_str(); 
}

//----------------------------------------
// Accelerators_C
//----------------------------------------
Accelerators_C::Accelerators_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions) 
	: Resource_C(pRCFile, pszName)
{
	m_sOptions = pszOptions;
}

Accelerators_C::~Accelerators_C()
{

}
	
const wchar_t* Accelerators_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void Accelerators_C::AddAccelerator(const wchar_t *pszAccelerator)
{
	m_accelerators.push_back(pszAccelerator);
}

//----------------------------------------
// Bitmap_C
//----------------------------------------
Bitmap_C::Bitmap_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName) 
	: Resource_C(pRCFile, pszName)
{
	m_sFileName = pszFileName;
}

Bitmap_C::~Bitmap_C()
{

}
	
const wchar_t* Bitmap_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

//----------------------------------------
// Cursor_C
//----------------------------------------
Cursor_C::Cursor_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName) 
	: Resource_C(pRCFile, pszName)
{
	m_sFileName = pszFileName;
}

Cursor_C::~Cursor_C()
{

}
	
const wchar_t* Cursor_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

//----------------------------------------
// Font_C
//----------------------------------------
Font_C::Font_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName) 
	: Resource_C(pRCFile, pszName)
{
	m_sFileName = pszFileName;
}

Font_C::~Font_C()
{

}
	
const wchar_t* Font_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

//----------------------------------------
// Icon_C
//----------------------------------------
Icon_C::Icon_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName) 
	: Resource_C(pRCFile, pszName)
{
	m_sFileName = pszFileName;
}

Icon_C::~Icon_C()
{

}
	
const wchar_t* Icon_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

//----------------------------------------
// RCData_C
//----------------------------------------
RCData_C::RCData_C(RCFile_C *pRCFile, const wchar_t *pszName) 
	: Resource_C(pRCFile, pszName)
{
}

RCData_C::~RCData_C()
{

}
	
const wchar_t* RCData_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

const wchar_t* RCData_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

const wchar_t* RCData_C::GetRawData() 
{ 
	return m_sData.c_str(); 
}

void RCData_C::SetFileName(const wchar_t *pszFileName)
{
	m_sFileName = pszFileName;
}

void RCData_C::SetOptions(const wchar_t *pszOptions)
{
	m_sOptions = pszOptions;
}

void RCData_C::SetRawData(const wchar_t *pszData)
{
	m_sData = pszData;
}

//----------------------------------------
// UserDef_C
//----------------------------------------
UserDef_C::UserDef_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszType) 
	: Resource_C(pRCFile, pszName)
{
	m_sType = pszType;
}

UserDef_C::~UserDef_C()
{

}
	
const wchar_t* UserDef_C::GetType() 
{ 
	return m_sType.c_str(); 
}

const wchar_t* UserDef_C::GetFileName() 
{ 
	return m_sFileName.c_str(); 
}

const wchar_t* UserDef_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void UserDef_C::SetFileName(const wchar_t *pszFileName)
{
	m_sFileName = pszFileName;
}

void UserDef_C::SetOptions(const wchar_t *pszOptions)
{
	m_sOptions = pszOptions;
}

//----------------------------------------
// Dialog_C
//----------------------------------------
Dialog_C::Dialog_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszCaption, size_t nCaptionStart, size_t nCaptionEnd) 
	: Resource_C(pRCFile, pszName)
{
	m_sCaption = pszCaption;
	m_nCaptionStartPos = nCaptionStart;
	m_nCaptionEndPos = nCaptionEnd;
	m_bTranslationHint = false;
}

Dialog_C::~Dialog_C()
{

}
	
const wchar_t* Dialog_C::GetCaption() 
{ 
	return m_sCaption.c_str(); 
}

void Dialog_C::AddControl(Control_C *pControl)
{
	m_controls.push_back(pControl);
}

int Dialog_C::GetControlsCount()
{
	return m_controls.size();
}

Control_C* Dialog_C::GetControl(int idx)
{
	if (idx < 0 || idx >= (int)m_controls.size())
		return NULL;
	else
		return m_controls[idx];
}

//----------------------------------------
// Control_C
//----------------------------------------
Control_C::Control_C(const wchar_t *pszType, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd)
{
	m_sType = pszType;
	m_sText = pszText;
	m_nTextStartPos = nTextStart;
	m_nTextEndPos = nTextEnd;
	m_sID = pszID;
	m_sOptions = pszOptions;
	m_bTranslationHint = false;
}

Control_C::~Control_C()
{

}
	
const wchar_t* Control_C::GetType() 
{ 
	return m_sType.c_str(); 
}

const wchar_t* Control_C::GetText() 
{ 
	return m_sText.c_str(); 
}

const wchar_t* Control_C::GetID() 
{ 
	return m_sID.c_str(); 
}

const wchar_t* Control_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

//----------------------------------------
// Menu_C
//----------------------------------------
Menu_C::Menu_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions) 
	: Resource_C(pRCFile, pszName)
{
	m_sOptions = pszOptions;
}

Menu_C::~Menu_C()
{

}
	
const wchar_t* Menu_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void Menu_C::AddItem(MenuItem_C *pItem)
{
	m_items.push_back(pItem);
}

int Menu_C::GetItemsCount()
{
	return m_items.size();
}

MenuItem_C* Menu_C::GetItem(int idx)
{
	if (idx < 0 || idx >= (int)m_items.size())
		return NULL;
	else
		return m_items[idx];
}

//----------------------------------------
// MenuItem_C
//----------------------------------------
MenuItem_C::MenuItem_C(const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd)
{
	m_sText = pszText;
	m_nTextStartPos = nTextStart;
	m_nTextEndPos = nTextEnd;
	m_sID = pszID;
	m_sOptions = pszOptions;
	m_bTranslationHint = false;
}

MenuItem_C::~MenuItem_C()
{

}
	
const wchar_t* MenuItem_C::GetText() 
{ 
	return m_sText.c_str(); 
}

const wchar_t* MenuItem_C::GetID() 
{ 
	return m_sID.c_str(); 
}

const wchar_t* MenuItem_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

//----------------------------------------
// MenuPopup_C
//----------------------------------------
MenuPopup_C::MenuPopup_C(MenuPopup_C *pOwner, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd)
	: MenuItem_C(pszText, pszID, pszOptions, nTextStart, nTextEnd)
{
	m_pOwner = pOwner;
}

MenuPopup_C::~MenuPopup_C()
{

}
	
MenuPopup_C* MenuPopup_C::GetOwner() 
{ 
	return m_pOwner; 
}

void MenuPopup_C::AddItem(MenuItem_C *pItem)
{
	m_items.push_back(pItem);
}

int MenuPopup_C::GetItemsCount()
{
	return m_items.size();
}

MenuItem_C* MenuPopup_C::GetItem(int idx)
{
	if (idx < 0 || idx >= (int)m_items.size())
		return NULL;
	else
		return m_items[idx];
}

//----------------------------------------
// StringTable_C
//----------------------------------------
StringTable_C::StringTable_C(RCFile_C *pRCFile, const wchar_t *pszOptions) 
	: Resource_C(pRCFile, L"")
{
	m_sOptions = pszOptions;
}

StringTable_C::~StringTable_C()
{

}
	
const wchar_t* StringTable_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void StringTable_C::AddString(String_C *pString)
{
	m_string.push_back(pString);
	m_pRCFile->AddString(pString->GetID(), pString);
}

int StringTable_C::GetStringsCount()
{
	return m_string.size();
}

String_C* StringTable_C::GetString(int idx)
{
	if (idx < 0 || idx >= (int)m_string.size())
		return NULL;
	else
		return m_string[idx];
}

//----------------------------------------
// String_C
//----------------------------------------
String_C::String_C(const wchar_t *pszID, const wchar_t *pszText, size_t nTextStart, size_t nTextEnd)
{
	m_sText = pszText;
	m_nTextStartPos = nTextStart;
	m_nTextEndPos = nTextEnd;
	m_sID = pszID;
	m_bTranslationHint = false;
}

String_C::~String_C()
{

}
	
const wchar_t* String_C::GetID() 
{ 
	return m_sID.c_str(); 
}

const wchar_t* String_C::GetText() 
{ 
	return m_sText.c_str(); 
}

//----------------------------------------
// VersionInfo_C
//----------------------------------------
VersionInfo_C::VersionInfo_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions) 
	: Resource_C(pRCFile, pszName)
{
	m_sOptions = pszOptions;
}

VersionInfo_C::~VersionInfo_C()
{

}
	
const wchar_t* VersionInfo_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void VersionInfo_C::AddBlock(VersionInfoBlock_C *pBlock)
{
	m_blocks.push_back(pBlock);
}


//----------------------------------------
// VersionInfoValue_C
//----------------------------------------
VersionInfoValue_C::VersionInfoValue_C(const wchar_t *pszID, const wchar_t *pszText)
{
	m_sText = pszText;
	m_sID = pszID;
}

VersionInfoValue_C::~VersionInfoValue_C()
{

}
	
const wchar_t* VersionInfoValue_C::GetID() 
{ 
	return m_sID.c_str(); 
}

const wchar_t* VersionInfoValue_C::GetText() 
{ 
	return m_sText.c_str(); 
}

//----------------------------------------
// VersionInfoBlock_C
//----------------------------------------
VersionInfoBlock_C::VersionInfoBlock_C(VersionInfoBlock_C *pOwner, const wchar_t *pszID)
{
	m_pOwner = pOwner;
	m_sID = pszID;
}

VersionInfoBlock_C::~VersionInfoBlock_C()
{

}
	
VersionInfoBlock_C* VersionInfoBlock_C::GetOwner() 
{ 
	return m_pOwner; 
}

const wchar_t* VersionInfoBlock_C::GetID() 
{ 
	return m_sID.c_str(); 
}

void VersionInfoBlock_C::AddBlock(VersionInfoBlock_C *pBlock)
{
	m_blocks.push_back(pBlock);
}

void VersionInfoBlock_C::AddValue(VersionInfoValue_C *pValue)
{
	m_values.push_back(pValue);
}

//----------------------------------------
// Toolbar_C
//----------------------------------------
Toolbar_C::Toolbar_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions) 
	: Resource_C(pRCFile, pszName)
{
	m_sOptions = pszOptions;
}

Toolbar_C::~Toolbar_C()
{

}
	
const wchar_t* Toolbar_C::GetOptions() 
{ 
	return m_sOptions.c_str(); 
}

void Toolbar_C::AddItem(const wchar_t *pszItem)
{
	m_items.push_back(pszItem);
}

int Toolbar_C::GetItemsCount()
{
	return m_items.size();
}

const wchar_t *Toolbar_C::GetItem(int idx)
{
	if (idx < 0 || idx >= (int)m_items.size())
		return NULL;
	else
		return m_items[idx].c_str();
}
