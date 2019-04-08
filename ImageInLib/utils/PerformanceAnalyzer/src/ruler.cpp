#include "stdafx.h"
#include "ruler.h"
#include "FormatTime.h"

#define MIN_SEGMENT_WIDTH_BASE 5.0
CBrush Ruler::m_rulerBrush(RGB(200, 200, 200));
CPen Ruler::m_rulerPen(PS_SOLID, 1, RGB(0,0,0));

Ruler::Ruler(CDC & dc, double minimum, double range, const CRect & rulerRect, int leftOffset) :
	m_dc(dc),
	m_min(minimum),
	m_range(range),
	m_rulerRect(rulerRect),
	m_leftOffset(leftOffset)
{
}

Ruler::~Ruler()
{
}

void Ruler::DrawSegments(bool bBig, double step, double leftPoint, double rightPoint, bool bInverse)
{
	int h = bBig ? 20 : 10;
	if (bInverse)
	{
		for (int idx = 0;; ++idx)
		{
			double posX = rightPoint - idx * step;
			if (posX < leftPoint)
				break;

			m_dc.MoveTo(int(posX + 0.9999999999), m_rulerRect.bottom - 2);
			m_dc.LineTo(int(posX + 0.9999999999), m_rulerRect.bottom - h);
		}
	}
	else
	{
		for (int idx = 0;; ++idx)
		{
			double posX = leftPoint + idx * step;
			if (posX > rightPoint)
				break;

			m_dc.MoveTo(int(posX + 0.99999999999), m_rulerRect.bottom - 2);
			m_dc.LineTo(int(posX + 0.99999999999), m_rulerRect.bottom - h);
		}
	}
}

void Ruler::Draw()
{
	double valueRange = m_range;
	
	CBrush * rulerBrushOld = m_dc.SelectObject(&m_rulerBrush);
	CPen * pOldPen = (CPen*)m_dc.SelectStockObject(NULL_PEN);
	CRect rect(m_rulerRect);

	// draw from beginning and eliminate null pen border
	rect.left -= m_leftOffset;
	rect.bottom += 1;
	rect.right += 1;

	m_dc.Rectangle(rect);
	m_dc.SelectObject(&rulerBrushOld);
	if (valueRange <= 0.0)
		return;

	m_dc.SelectObject(&m_rulerPen);

	double rulerBaseWidth = m_rulerRect.Width() - 1;
	int segmentCountMax = int(rulerBaseWidth / MIN_SEGMENT_WIDTH_BASE);
	double unit = 100000;
	int segmentCount = 0;
	double segmentOffsetInPx, segmentSizeInPx;
	segmentOffsetInPx = segmentSizeInPx = 0.0;
	while (true)
	{
		if (unit > (valueRange / double(segmentCountMax)))
		{
			unit /= 10.0;
		}
		else
		{
			while (true)
			{
				unit *= 10.0;
				segmentCount = (int)((m_range / unit) + 0.5);
				double offset = ((trunc(m_min / unit) + 1) * unit) - m_min;
				segmentOffsetInPx = m_rulerRect.left + ((double)rulerBaseWidth * offset / valueRange);
				segmentSizeInPx = ((double)rulerBaseWidth * unit) / valueRange;

				if (abs(segmentSizeInPx - (segmentOffsetInPx - m_rulerRect.left)) < 0.000001)
				{

					segmentOffsetInPx = m_rulerRect.left;
					segmentCount++;
				}

				if (((double)segmentSizeInPx / 4.0) < MIN_SEGMENT_WIDTH_BASE)
					continue;

				break;
			}

			break;
		}
	}

	double smallSegmentSize = segmentSizeInPx / 4.0;
	DrawSegments(false, smallSegmentSize, m_rulerRect.left, segmentOffsetInPx, true);
	DrawSegments(true, segmentSizeInPx, segmentOffsetInPx, segmentOffsetInPx + segmentCount * segmentSizeInPx, false);
	for (int i = 0; i < segmentCount; i++)
	{
		DrawSegments(false, smallSegmentSize, segmentOffsetInPx, segmentOffsetInPx + segmentSizeInPx, false);
		segmentOffsetInPx += segmentSizeInPx;
	}

	m_dc.SelectObject(pOldPen);

	CString units = L"x" + FormatTimeInterval(unit, 0);
	CSize sz = m_dc.GetTextExtent(units);
	m_dc.TextOut(m_rulerRect.right - sz.cx - 5, m_rulerRect.top + 5, units);
}