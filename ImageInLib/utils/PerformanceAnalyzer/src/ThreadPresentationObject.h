#pragma once

#include <vector>
#include <map>
#include "FunctionPresentationObject.h"

class CThreadPresentationObject
{
public:
	CThreadPresentationObject(LONG yPosition, LONG nHeaderWidth, COLORREF clrHeader, LPCTSTR lpszHeaderText, COLORREF clrBackground);
	~CThreadPresentationObject();

	void AddFunctions(std::map<CString, COLORREF> & clrFunctions, const std::vector<const DataModel::CFunction*>& functions);
	void Invalidate();

	int Draw(CDC &dc, CSize size, CPoint ptPositionShift, double xScrollPos, int yScrollPos);

	LONG Height() const { return m_szHeader.cy; }

	bool GetTooltipText(CPoint point, CString& strTooltip, int yScrollPos) const;

protected:
	int      m_yPosition;
	CSize    m_szHeader;
	CString  m_strHeaderText;
	COLORREF m_clrHeader;
	COLORREF m_clrBackground;

	double   m_dTimeOffset;
	double   m_dMinTime;
	double   m_dMaxTime;

	// list of function presentation objects sorted by start position
	std::vector<CFunctionPresentationObject> m_functions;
};

