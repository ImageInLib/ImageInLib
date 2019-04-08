// rcfile.h: interface for the RCFile_C class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_)
#define AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class RCFile_C;

class Resource_C 
{
public:
	Resource_C(RCFile_C *pRCFile, const wchar_t *pszName);
	virtual ~Resource_C();
	const wchar_t *GetName();
protected:
	RCFile_C *m_pRCFile;
	std::wstring m_sName;
};


class Accelerators_C : public Resource_C
{
public:
	Accelerators_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~Accelerators_C();
	const wchar_t *GetOptions();
	void AddAccelerator(const wchar_t *pszAccelerator);
protected:
	std::wstring m_sOptions;
	std::vector<std::wstring> m_accelerators;
};


class Bitmap_C : public Resource_C
{
public:
	Bitmap_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Bitmap_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Cursor_C : public Resource_C
{
public:
	Cursor_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Cursor_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Font_C : public Resource_C
{
public:
	Font_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Font_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Icon_C : public Resource_C
{
public:
	Icon_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Icon_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class RCData_C : public Resource_C
{
public:
	RCData_C(RCFile_C *pRCFile, const wchar_t *pszName); 
	virtual ~RCData_C();
	const wchar_t *GetFileName();
	const wchar_t *GetOptions();
	const wchar_t *GetRawData();
	void SetFileName(const wchar_t *pszFileName);
	void SetOptions(const wchar_t *pszOptions);
	void SetRawData(const wchar_t *pszData);
protected:
	std::wstring m_sFileName;
	std::wstring m_sOptions;
	std::wstring m_sData;
};


class Control_C
{
public:
	Control_C(const wchar_t *pszType, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd); 
	virtual ~Control_C();
	const wchar_t *GetType();
	const wchar_t *GetText();
	size_t GetTextStartPos() { return m_nTextStartPos; }
	size_t GetTextEndPos() { return m_nTextEndPos; }
	const wchar_t *GetID();
	const wchar_t *GetOptions();
	const wchar_t *GetTranslationHint() { return m_bTranslationHint ? m_sTranslationHint.c_str() : NULL; }
	void SetTranslationHint(const wchar_t* pszHint) { m_bTranslationHint = pszHint ? true : false; if (m_bTranslationHint) m_sTranslationHint = pszHint; }
protected:
	std::wstring m_sType;
	size_t m_nTextStartPos;
	size_t m_nTextEndPos;
	std::wstring m_sText;
	bool m_bTranslationHint;
	std::wstring m_sTranslationHint;
	std::wstring m_sID;
	std::wstring m_sOptions;
};


class Dialog_C : public Resource_C
{
public:
	Dialog_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszCaption, size_t nCaptionStart, size_t nCaptionEnd);
	virtual ~Dialog_C();
	const wchar_t *GetCaption();
	size_t GetCaptionStartPos() { return m_nCaptionStartPos; }
	size_t GetCaptionEndPos() { return m_nCaptionEndPos; }
	void AddControl(Control_C *pControl);
	int GetControlsCount();
	Control_C *GetControl(int idx);
	const wchar_t *GetTranslationHint() { return m_bTranslationHint ? m_sTranslationHint.c_str() : NULL; }
	void SetTranslationHint(const wchar_t* pszHint) { m_bTranslationHint = pszHint ? true : false; if (m_bTranslationHint) m_sTranslationHint = pszHint; }
protected:
	size_t m_nCaptionStartPos;
	size_t m_nCaptionEndPos;
	std::wstring m_sCaption;
	bool m_bTranslationHint;
	std::wstring m_sTranslationHint;
	std::vector<Control_C*> m_controls;
};


class MenuItem_C
{
public:
	MenuItem_C(const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd); 
	virtual ~MenuItem_C();
	const wchar_t *GetText();
	size_t GetTextStartPos() { return m_nTextStartPos; }
	size_t GetTextEndPos() { return m_nTextEndPos; }
	const wchar_t *GetID();
	const wchar_t *GetOptions();
	const wchar_t *GetTranslationHint() { return m_bTranslationHint ? m_sTranslationHint.c_str() : NULL; }
	void SetTranslationHint(const wchar_t* pszHint) { m_bTranslationHint = pszHint ? true : false; if (m_bTranslationHint) m_sTranslationHint = pszHint; }
protected:
	std::wstring m_sText;
	size_t m_nTextStartPos;
	size_t m_nTextEndPos;
	bool m_bTranslationHint;
	std::wstring m_sTranslationHint;
	std::wstring m_sID;
	std::wstring m_sOptions;
};


class MenuPopup_C : public MenuItem_C
{
public:
	MenuPopup_C(MenuPopup_C *pOwner, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions, size_t nTextStart, size_t nTextEnd); 
	virtual ~MenuPopup_C();
	MenuPopup_C *GetOwner();
	void AddItem(MenuItem_C *pControl);
	int GetItemsCount();
	MenuItem_C *GetItem(int idx);
protected:
	MenuPopup_C *m_pOwner;
	std::vector<MenuItem_C*> m_items;
};


class Menu_C : public Resource_C
{
public:
	Menu_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~Menu_C();
	const wchar_t *GetOptions();
	void AddItem(MenuItem_C *pControl);
	int GetItemsCount();
	MenuItem_C *GetItem(int idx);
protected:
	std::wstring m_sOptions;
	std::vector<MenuItem_C*> m_items;
};


class String_C
{
public:
	String_C(const wchar_t *pszID, const wchar_t *pszText, size_t nTextStart, size_t nTextEnd); 
	virtual ~String_C();
	const wchar_t *GetID();
	const wchar_t *GetText();
	size_t GetTextStartPos() { return m_nTextStartPos; }
	size_t GetTextEndPos() { return m_nTextEndPos; }
	const wchar_t *GetTranslationHint() { return m_bTranslationHint ? m_sTranslationHint.c_str() : NULL; }
	void SetTranslationHint(const wchar_t* pszHint) { m_bTranslationHint = pszHint ? true : false; if (m_bTranslationHint) m_sTranslationHint = pszHint; }
protected:
	std::wstring m_sID;
	std::wstring m_sText;
	size_t m_nTextStartPos;
	size_t m_nTextEndPos;
	bool m_bTranslationHint;
	std::wstring m_sTranslationHint;
};


class StringTable_C : public Resource_C
{
public:
	StringTable_C(RCFile_C *pRCFile, const wchar_t *pszOptions); 
	virtual ~StringTable_C();
	const wchar_t *GetOptions();
	void AddString(String_C *pString);
	int GetStringsCount();
	String_C *GetString(int idx);
protected:
	std::wstring m_sOptions;
	std::vector<String_C*> m_string;
};


class VersionInfoValue_C
{
public:
	VersionInfoValue_C(const wchar_t *pszID, const wchar_t *pszText); 
	virtual ~VersionInfoValue_C();
	const wchar_t *GetID();
	const wchar_t *GetText();
protected:
	std::wstring m_sID;
	std::wstring m_sText;
};


class VersionInfoBlock_C
{
public:
	VersionInfoBlock_C(VersionInfoBlock_C *pOwner, const wchar_t *pszID); 
	virtual ~VersionInfoBlock_C();
	const wchar_t *GetID();
	VersionInfoBlock_C *GetOwner();
	void AddBlock(VersionInfoBlock_C *pBlock);
	void AddValue(VersionInfoValue_C *pValue);
protected:
	VersionInfoBlock_C *m_pOwner;
	std::wstring m_sID;
	std::vector<VersionInfoBlock_C*> m_blocks;
	std::vector<VersionInfoValue_C*> m_values;
};


class VersionInfo_C : public Resource_C
{
public:
	VersionInfo_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~VersionInfo_C();
	const wchar_t *GetOptions();
	void AddBlock(VersionInfoBlock_C *pBlock);
protected:
	std::wstring m_sOptions;
	std::vector<VersionInfoBlock_C*> m_blocks;
};


class Toolbar_C : public Resource_C
{
public:
	Toolbar_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~Toolbar_C();
	const wchar_t *GetOptions();
	void AddItem(const wchar_t *pszItem);
	int GetItemsCount();
	const wchar_t *GetItem(int idx);
protected:
	std::wstring m_sOptions;
	std::vector<std::wstring> m_items;
};


class UserDef_C : public Resource_C
{
public:
	UserDef_C(RCFile_C *pRCFile, const wchar_t *pszName, const wchar_t *pszType); 
	virtual ~UserDef_C();
	const wchar_t *GetFileName();
	const wchar_t *GetOptions();
	const wchar_t *GetType();
	void SetFileName(const wchar_t *pszFileName);
	void SetOptions(const wchar_t *pszOptions);
protected:
	std::wstring m_sType;
	std::wstring m_sFileName;
	std::wstring m_sOptions;
};


class RCFile_C  
{
public:
	RCFile_C();
	virtual ~RCFile_C();

	void AddResource(Resource_C *pResource);

	int GetResourceCount();
	Resource_C *GetResource(int idx);
	int GetResource(int idx, const type_info *type, Resource_C **ppResource);

	void AddString(const wchar_t* stringID, String_C *pString);

	String_C* FindString(const wchar_t* stringID);
	Menu_C* FindMenu(const wchar_t* menuID);
	Dialog_C* FindDialog(const wchar_t* dialogID);

protected:
	std::vector<Resource_C*> m_resources;
	std::map<std::wstring, String_C*> m_strings;
	std::map<std::wstring, Menu_C*> m_menus;
	std::map<std::wstring, Dialog_C*> m_dialogs;
};

#endif // !defined(AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_)
