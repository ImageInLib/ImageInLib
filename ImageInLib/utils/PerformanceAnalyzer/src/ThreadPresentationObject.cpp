#include "stdafx.h"
#include "ThreadPresentationObject.h"
#include "ThreadContainer.h"
#include "ThreadsList.h"
#include "Global.h"


CThreadPresentationObject::CThreadPresentationObject(
	LONG yPosition,
	LONG nHeaderWidth,
	COLORREF clrHeader,
	LPCTSTR lpszHeaderText,
	COLORREF clrBackground)
	:
	m_yPosition(yPosition),
	m_strHeaderText(lpszHeaderText),
	m_clrBackground(clrBackground)
{
	// set header size and color
	m_szHeader.cx = nHeaderWidth;
	m_szHeader.cy = CFunctionPresentationObject::LevelHeight() + CFunctionPresentationObject::LevelSpacing();
	m_clrHeader = clrHeader;
}


CThreadPresentationObject::~CThreadPresentationObject()
{
}

void CThreadPresentationObject::AddFunctions(std::map<CString, COLORREF> & clrFunctions, const std::vector<const DataModel::CFunction*>& functions)
{
	UINT nLevel = 0;
	m_functions.reserve(m_functions.size() + functions.size());
	for (auto function : functions)
	{
		CFunctionPresentationObject newFunction(function, clrFunctions);
		m_functions.push_back(newFunction);
		if (function->Level > nLevel)
			nLevel = function->Level;
	}
	// update header height
	LONG nHight = (nLevel + 1) * (CFunctionPresentationObject::LevelHeight() + CFunctionPresentationObject::LevelSpacing());
	if (nHight > m_szHeader.cy)
		m_szHeader.cy = nHight;
}

void CThreadPresentationObject::Invalidate()
{
	// update size of each function presentation object
	for (CFunctionPresentationObject& function : m_functions)
	{
		function.Invalidate();
	}
}

bool CThreadPresentationObject::GetTooltipText(CPoint point, CString& strTooltip, int yScrollPos) const
{
	point.y += yScrollPos;
	// check if point is within thread object
	if (point.x < 0 || point.y < m_yPosition || point.y > m_yPosition + Height() - 1)
		return false;
	// check if point is within header rectangle
	if (point.x >= 0 && point.x < m_szHeader.cx) {
		if (m_strHeaderText.IsEmpty())
			strTooltip.Empty();
		else
			strTooltip = _T("Thread: ") + m_strHeaderText; // or just strTooltip.Empty();
		return true;
	}
	// check function objects
	point.y -= m_yPosition;
	point.x -= m_szHeader.cx;
	for (const CFunctionPresentationObject& function : m_functions)
	{
		if (function.GetTooltipText(point, strTooltip))
			return true;
	}
	return true;
}

int CThreadPresentationObject::Draw(CDC &dc, CSize size, CPoint ptPositionShift, double xScrollPos, int yScrollPos)
{
	ptPositionShift.y -= yScrollPos;
	// check if thread is visible in the area specified by the rectArea
	// if we are not in the area yet, skip the thread
	if (m_yPosition + Height() <= yScrollPos)
		return -1;
	// if we are out of range, end the loop
	if (m_yPosition > yScrollPos + size.cy)
		return 1;

	CRect rectTPO(CPoint(0, m_yPosition), CSize(size.cx, Height()));

	CBrush *pOldBrush;
	CPen   *pOldPen;

	// fill background area for under function presentation objects
	CRect rectBkgnd(rectTPO);
	rectBkgnd.left += m_szHeader.cx;
	rectBkgnd += ptPositionShift;
	CPen bkgndPen(PS_SOLID, 1, m_clrBackground);
	CBrush bkgndBrush(m_clrBackground);
	pOldBrush = dc.SelectObject(&bkgndBrush);
	pOldPen = dc.SelectObject(&bkgndPen);
	dc.Rectangle(rectBkgnd);

	// draw the header (if present)
	if (m_szHeader.cx > 0) {
		CRect rectHeader(rectTPO);
		rectHeader.right = rectHeader.left + m_szHeader.cx;
		rectHeader += ptPositionShift;
		CPen hdrPen(PS_SOLID, 1, m_clrHeader);
		CBrush hdrBrush(m_clrHeader);
		dc.SelectObject(&hdrBrush);
		dc.SelectObject(&hdrPen);
		dc.Rectangle(rectHeader);
		dc.DrawText(m_strHeaderText, rectHeader, DT_CENTER | DT_VCENTER | DT_SINGLELINE);
	}

	CRgn rgn;
	rgn.CreateRectRgnIndirect(rectBkgnd);
	dc.SelectClipRgn(&rgn, RGN_AND);

	// draw funciton presentation objects
	ptPositionShift.x += m_szHeader.cx;
	ptPositionShift.y += m_yPosition;
	for (auto& function : m_functions)
	{
		if (function.Draw(dc, xScrollPos, ptPositionShift, rectBkgnd.Width()) > 0)
			break;
	}

	dc.SelectObject(pOldPen);
	dc.SelectObject(pOldBrush);

	return 0;
}
