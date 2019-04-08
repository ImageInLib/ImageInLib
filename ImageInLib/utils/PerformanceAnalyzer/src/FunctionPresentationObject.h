#pragma once

#include "Function.h"
#include "ZoomInfo.h"
#include <map>

class CFunctionPresentationObject
{
public:
	CFunctionPresentationObject(const DataModel::CFunction* pFunctionData, std::map<CString, COLORREF> & clrFunctions);
	~CFunctionPresentationObject();

	static int LevelHeight() { return 30; }
	static int LevelSpacing() { return 1; }

	void Invalidate();
	void SetColor(COLORREF color);
	bool GetTooltipText(CPoint point, CString& strTooltip) const;

	int Draw(CDC &dc, double dScrollPos, CPoint ptPositionShift, int nDisplayingAreaWidth);

	static void ExtractFunctionName(const CString& strSrc, CString& strDst);

	static COLORREF GetNextColor_FunctionPresentationObject();
	static void InitColorIndexes();

private:
	void UpdateTextForDisplaying(CDC& dc, double dObjectWidth);

	int Top() const { return m_pFunctionData->Level * (LevelHeight() + LevelSpacing()); }
	int Bottom() const { return Top() + LevelHeight(); }

protected:
	CRect    m_rect;
	COLORREF m_clrCore;
	COLORREF m_clrBorder;
	bool     m_bTextIsValid;
	CString  m_text;
	CRect    m_rectText;
	// reference to data
	const DataModel::CFunction* m_pFunctionData;
	static size_t m_iColorFPO;
};
