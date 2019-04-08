#include "TimeMark.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define CURSOR_COLOR  RGB(255,0,0)

#define TIME_NOT_SET -1.0

TimeMark::TimeMark(DataModel::CZoomInfo& zoomInfo)
	: m_zoomInfo(zoomInfo)
{
	m_dTime = TIME_NOT_SET;
}

void TimeMark::Draw(CDC* pDC, CRect rect)
{
	if (!pDC || m_dTime == TIME_NOT_SET )
		return;

	CPen cursorPen(PS_SOLID, 1, CURSOR_COLOR);
	CGdiObject* pPenOld = pDC->SelectObject(&cursorPen);

	CPoint ptPos = rect.TopLeft();
	double dTimeOffset = max(0, m_dTime - m_zoomInfo.ZoomRangeStart() + m_zoomInfo.DataShift());
	ptPos.x += (int)(m_zoomInfo.TimePositionToPixelPosition(dTimeOffset)+0.5);
	pDC->MoveTo(ptPos);
	ptPos.y = rect.bottom;
	pDC->LineTo(ptPos);

	pDC->SelectObject(pPenOld);
}

bool TimeMark::SetTime(double dTime)
{
	if (dTime < 0.0)
		dTime = 0.0;

	if (dTime > m_zoomInfo.DataWidth())
		dTime = m_zoomInfo.DataWidth();

	if (m_dTime == dTime)
		return false;

	m_dTime = dTime;
	return true;
}