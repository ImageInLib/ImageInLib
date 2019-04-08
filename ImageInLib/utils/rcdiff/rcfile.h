// rcfile.h: interface for the RCFile_C class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_)
#define AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class Resource_C 
{
public:
	Resource_C(const wchar_t *pszName);
	virtual ~Resource_C();
	const wchar_t *GetName();
protected:
	std::wstring m_sName;
};


class Accelerators_C : public Resource_C
{
public:
	Accelerators_C(const wchar_t *pszName, const wchar_t *pszOptions); 
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
	Bitmap_C(const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Bitmap_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Cursor_C : public Resource_C
{
public:
	Cursor_C(const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Cursor_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Font_C : public Resource_C
{
public:
	Font_C(const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Font_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class Icon_C : public Resource_C
{
public:
	Icon_C(const wchar_t *pszName, const wchar_t *pszFileName); 
	virtual ~Icon_C();
	const wchar_t *GetFileName();
protected:
	std::wstring m_sFileName;
};


class RCData_C : public Resource_C
{
public:
	RCData_C(const wchar_t *pszName); 
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
	Control_C(const wchar_t *pszType, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions); 
	virtual ~Control_C();
	const wchar_t *GetType();
	const wchar_t *GetText();
	const wchar_t *GetID();
	const wchar_t *GetOptions();
protected:
	std::wstring m_sType;
	std::wstring m_sText;
	std::wstring m_sID;
	std::wstring m_sOptions;
};


class Dialog_C : public Resource_C
{
public:
	Dialog_C(const wchar_t *pszName, const wchar_t *pszType, const wchar_t *pszOptions); 
	virtual ~Dialog_C();
	const wchar_t *GetOptions();
	const wchar_t *GetType();
	void AddControl(Control_C *pControl);
	int GetControlsCount();
	Control_C *GetControl(int idx);
protected:
	std::wstring m_sType;
	std::wstring m_sOptions;
	vector<Control_C*> m_controls;
};


class MenuItem_C
{
public:
	MenuItem_C(const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions); 
	virtual ~MenuItem_C();
	const wchar_t *GetText();
	const wchar_t *GetID();
	const wchar_t *GetOptions();
protected:
	std::wstring m_sText;
	std::wstring m_sID;
	std::wstring m_sOptions;
};


class MenuPopup_C : public MenuItem_C
{
public:
	MenuPopup_C(MenuPopup_C *pOwner, const wchar_t *pszText, const wchar_t *pszID, const wchar_t *pszOptions); 
	virtual ~MenuPopup_C();
	MenuPopup_C *GetOwner();
	const wchar_t *GetText();
	const wchar_t *GetID();
	const wchar_t *GetOptions();
	void AddItem(MenuItem_C *pControl);
	int GetItemsCount();
	MenuItem_C *GetItem(int idx);
protected:
	MenuPopup_C *m_pOwner;
	vector<MenuItem_C*> m_items;
};


class Menu_C : public Resource_C
{
public:
	Menu_C(const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~Menu_C();
	const wchar_t *GetOptions();
	void AddItem(MenuItem_C *pControl);
	int GetItemsCount();
	MenuItem_C *GetItem(int idx);
protected:
	std::wstring m_sOptions;
	vector<MenuItem_C*> m_items;
};


class String_C
{
public:
	String_C(const wchar_t *pszID, const wchar_t *pszText); 
	virtual ~String_C();
	const wchar_t *GetID();
	const wchar_t *GetText();
protected:
	std::wstring m_sID;
	std::wstring m_sText;
};


class StringTable_C : public Resource_C
{
public:
	StringTable_C(const wchar_t *pszOptions); 
	virtual ~StringTable_C();
	const wchar_t *GetOptions();
	void AddString(String_C *pString);
	int GetStringsCount();
	String_C *GetString(int idx);
protected:
	std::wstring m_sOptions;
	vector<String_C*> m_string;
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
	vector<VersionInfoBlock_C*> m_blocks;
	vector<VersionInfoValue_C*> m_values;
};


class VersionInfo_C : public Resource_C
{
public:
	VersionInfo_C(const wchar_t *pszName, const wchar_t *pszOptions); 
	virtual ~VersionInfo_C();
	const wchar_t *GetOptions();
	void AddBlock(VersionInfoBlock_C *pBlock);
protected:
	std::wstring m_sOptions;
	vector<VersionInfoBlock_C*> m_blocks;
};


class Toolbar_C : public Resource_C
{
public:
	Toolbar_C(const wchar_t *pszName, const wchar_t *pszOptions); 
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
	UserDef_C(const wchar_t *pszName, const wchar_t *pszType); 
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

protected:
	vector<Resource_C*> m_resources;
};

#endif // !defined(AFX_RCFILE_H__9D3078B8_6C5A_4913_89B9_BF47B04B0BB0__INCLUDED_)
