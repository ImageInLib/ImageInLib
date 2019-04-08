#include "TimeMark.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


void TimeMark::Draw(CDC* pDC, CRect rect, Style style)
{
	if (!pDC || m_dTime == TIME_NOT_SET )
		return;

	CPoint ptPos = rect.TopLeft();
	double dTimeOffset = m_dTime - m_zoomInfo.ZoomRangeStart() + m_zoomInfo.DataShift();
	ptPos.x += (int)(m_zoomInfo.TimePositionToPixelPosition(dTimeOffset)+0.5);

	if (ptPos.x >= rect.left && ptPos.x < rect.right)
	{
		CPen cursorPen(PS_SOLID, 1, m_color);
		CGdiObject* pPenOld = pDC->SelectObject(&cursorPen);

		pDC->MoveTo(ptPos);
		ptPos.y = rect.bottom;
		pDC->LineTo(ptPos);

		if (style == LeftArrow || style == RightArrow)
		{
			CPoint arrow[3];
			arrow[0] = CPoint(ptPos.x, rect.top);
			arrow[1] = CPoint(ptPos.x, rect.top + 5);
			arrow[2] = CPoint(ptPos.x + (style == LeftArrow ? -5 : 5), rect.top);
			CBrush arrowBrush(m_color);
			CGdiObject* pBrushOld = pDC->SelectObject(&arrowBrush);
			pDC->Polygon(arrow, _countof(arrow));
			pDC->SelectObject(pBrushOld);
		}

		pDC->SelectObject(pPenOld);
	}
}

bool TimeMark::SetTime(double dTime)
{
	if (dTime < 0.0 && dTime != TIME_NOT_SET)
		dTime = 0.0;

	if (dTime > m_zoomInfo.DataWidth())
		dTime = m_zoomInfo.DataWidth();

	if (m_dTime == dTime)
		return false;

	m_dTime = dTime;
	return true;
}