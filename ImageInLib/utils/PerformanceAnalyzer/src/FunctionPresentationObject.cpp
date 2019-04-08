#include "stdafx.h"
#include "FunctionPresentationObject.h"
#include "FormatTime.h"
#include "Global.h"

static int MINIMAL_FUNCTION_OBJ_WIDTH = 1;

size_t CFunctionPresentationObject::m_iColorFPO = 0;

CFunctionPresentationObject::CFunctionPresentationObject(const DataModel::CFunction* pFunctionData, std::map<CString, COLORREF> & clrFunctions) :
	m_pFunctionData(pFunctionData),
	m_bTextIsValid(false)
{
	COLORREF clrFunction;
	auto itColor = clrFunctions.find(pFunctionData->Description);
	if (itColor == clrFunctions.end())
	{
		clrFunction = GetNextColor_FunctionPresentationObject();
		clrFunctions[pFunctionData->Description] = clrFunction;
}
	else
{
		clrFunction = itColor->second;
}
	SetColor(clrFunction);
}

CFunctionPresentationObject::~CFunctionPresentationObject()
{
}

void CFunctionPresentationObject::Invalidate()
{
	// invalidate rectangle's size
	m_rect.SetRectEmpty();
	// invalidate text
	m_bTextIsValid = false;
}

void CFunctionPresentationObject::SetColor(COLORREF color)
{
	m_clrCore = color;
	m_clrBorder = RGB(GetRValue(m_clrCore) / 2, GetGValue(m_clrCore) / 2, GetBValue(m_clrCore) / 2);
}

bool CFunctionPresentationObject::GetTooltipText(CPoint point, CString& strTooltip) const
{
	if (m_rect.PtInRect(point)) {
		strTooltip.Format(_T("Thread: %d<br>%s\nStart time: %s<br>Duration: %s"),
			m_pFunctionData->ThreadID,
			(LPCTSTR)m_pFunctionData->Description,
			FormatAbsoluteTime(m_pFunctionData->Start, DefaultPrecisionAfterSeconds),
			FormatTimeInterval(m_pFunctionData->Stop - m_pFunctionData->Start, DefaultIntervalPrecision)
			);
		return true;
	}
	return false;
}

void CFunctionPresentationObject::ExtractFunctionName(const CString& strSrc, CString& strDst)
{
	int iStart = strSrc.Find(_T("::"));
	if (iStart != -1) {
		iStart += 2;	// add ::
		strDst = strSrc.Right(strSrc.GetLength() - iStart);
	}
	else
		strDst.Empty();
}

void CFunctionPresentationObject::UpdateTextForDisplaying(CDC& dc, double dObjectWidth)
{
	dObjectWidth += 0.5;
	int nObjectWidth = (dObjectWidth > INT_MAX) ? INT_MAX : (int)dObjectWidth;

	m_rectText.SetRectEmpty();
	m_text.Empty();
	if (!m_pFunctionData->Description.IsEmpty() && nObjectWidth > 40) {
		m_text = m_pFunctionData->Description;
		dc.DrawText(m_text, m_rectText, DT_CALCRECT | DT_SINGLELINE);
		if (m_rectText.Width() >= nObjectWidth - 2) {
			ExtractFunctionName(m_pFunctionData->Description, m_text);
			m_rectText.SetRectEmpty();
			dc.DrawText(m_text, m_rectText, DT_CALCRECT | DT_SINGLELINE);
			if (m_rectText.Width() >= nObjectWidth - 2) {
				m_rectText.SetRect(0, 0, nObjectWidth - 12, LevelHeight() - 4);
				dc.DrawTextEx(m_text.GetBuffer(m_text.GetLength()), m_text.GetLength(), m_rectText, DT_CALCRECT | DT_SINGLELINE | DT_END_ELLIPSIS | DT_MODIFYSTRING, NULL);
				m_text.ReleaseBuffer();
			}
		}
	}
	m_bTextIsValid = true;
}

int CFunctionPresentationObject::Draw(CDC &dc, double dScrollPos, CPoint ptPositionShift, int nDisplayingAreaWidth)
{
	DataModel::CZoomInfo & zoomInfo = Global_C::Instance()->ZoomInfo();
	dScrollPos += zoomInfo.DataShift();
	double dRangeEnd = zoomInfo.ZoomRange() + dScrollPos;

	// check if function object is visible (is within displaying range)
	// if we are not in the range yet, skip the function
	if (m_pFunctionData->Stop < dScrollPos)
		return -1;
	// if we are already out of range, end the loop
	if (m_pFunctionData->Start > dRangeEnd)
		return 1;


	// transform to pixel values
	double dTmp = zoomInfo.TimePositionToPixelPosition(m_pFunctionData->Start - dScrollPos);
	if (dTmp < 0.)
		dTmp = -1;
	else
		dTmp += 0.5;
	m_rect.left = (int)dTmp;
	dTmp = trunc(zoomInfo.TimePositionToPixelPosition(m_pFunctionData->Stop - dScrollPos) + 0.5);
	if (dTmp > (double)nDisplayingAreaWidth)
		dTmp = nDisplayingAreaWidth + 1.;
	m_rect.right = (int)dTmp;
	m_rect.top =  Top();
	m_rect.bottom = m_rect.top + LevelHeight();

	if (m_rect.Width() < MINIMAL_FUNCTION_OBJ_WIDTH) {
		m_rect.right = m_rect.left + MINIMAL_FUNCTION_OBJ_WIDTH;
		if (zoomInfo.PixelPositionToTimePosition(m_rect.right) > zoomInfo.DataWidth())
			m_rect.OffsetRect(-1, 0);
	}

	// update text that will be displayed if it is not valid
	if (!m_bTextIsValid) {
		double dObjectWidth = zoomInfo.TimePositionToPixelPosition(m_pFunctionData->Stop - m_pFunctionData->Start);
		UpdateTextForDisplaying(dc, dObjectWidth);
	}

	// calc rectangle for text
	CRect rectText;
	if (!m_text.IsEmpty())
	{
		rectText.top = m_rect.top + m_rect.Height() / 2 - m_rectText.Height() / 2 - 1;
		rectText.bottom = rectText.top + m_rectText.Height();
		rectText.left = max(0, m_rect.left);
		rectText.right = min(nDisplayingAreaWidth, m_rect.right);
		if (rectText.Width() < m_rectText.Width())
			rectText.left = (m_rect.left < 0) ? rectText.right - m_rectText.Width() - 1 : m_rect.left + 1;
		else
			rectText.left = rectText.left + rectText.Width() / 2 - m_rectText.Width() / 2;
		rectText.right = rectText.left + m_rectText.Width();
	}

	// add position shift
	CRect rectObject(m_rect);
	rectObject.OffsetRect(ptPositionShift);
	rectText.OffsetRect(ptPositionShift);

	// draw the object
	CPen  pen(PS_SOLID, 1, m_clrBorder);
	CPen* pOldPen = dc.SelectObject(&pen);
	CBrush  brush(m_clrCore);
	CBrush* pOldBrush = dc.SelectObject(&brush);
	dc.Rectangle(rectObject);	// alternatively you can use: dc.RoundRect(m_rect, CPoint(5, 5));
	if (!m_text.IsEmpty()) {
		COLORREF textColor, oldTextColor;
		if (GetRValue(m_clrCore) + GetGValue(m_clrCore) + GetBValue(m_clrCore) > 3 * 128)
			textColor = RGB(0, 0, 0);
		else
			textColor = RGB(255, 255, 255);
		oldTextColor = dc.SetTextColor(textColor);
		dc.DrawText(m_text, rectText, DT_SINGLELINE);
		dc.SetTextColor(oldTextColor);
	}
	dc.SelectObject(pOldBrush);
	dc.SelectObject(pOldPen);
	return 0;
}

COLORREF CFunctionPresentationObject::GetNextColor_FunctionPresentationObject()
{
#if 0
	const COLORREF colors[] =
	{
		RGB(0xf0, 0xf8, 0xff),
		RGB(0xfa, 0x80, 0x72),
		RGB(0x87, 0xce, 0xfa),
		RGB(0xfa, 0xeb, 0xd7),
		RGB(0x7f, 0xff, 0xd4),
		RGB(0xf5, 0xf5, 0xdc),
		RGB(0xe6, 0xe6, 0xfa),
		RGB(0x9a, 0xcd, 0x32)
	};
	m_iColorFPO = (m_iColorFPO + 1) % _countof(colors);
	return colors[m_iColorFPO];
#else
	const COLORREF colors[] =
	{
		RGB(0xFF, 0x00, 0x00),
		RGB(0xFF, 0x40, 0x00),
		RGB(0xFF, 0x80, 0x00),
		RGB(0xFF, 0xC0, 0x00),
		RGB(0xFF, 0xFF, 0x00),
		RGB(0xFF, 0x00, 0x40),
		RGB(0xFF, 0x40, 0x40),
		RGB(0xFF, 0x80, 0x40),
		RGB(0xFF, 0xC0, 0x40),
		RGB(0xFF, 0xFF, 0x40),
		RGB(0xFF, 0x00, 0x80),
		RGB(0xFF, 0x40, 0x80),
		RGB(0xFF, 0x80, 0x80),
		RGB(0xFF, 0xC0, 0x80),
		RGB(0xFF, 0xFF, 0x80),
		RGB(0xFF, 0x00, 0xC0),
		RGB(0xFF, 0x40, 0xC0),
		RGB(0xFF, 0x80, 0xC0),
		RGB(0xFF, 0xC0, 0xC0),
		RGB(0xFF, 0xFF, 0xC0),
		RGB(0xFF, 0x00, 0xFF),
		RGB(0xFF, 0x40, 0xFF),
		RGB(0xFF, 0x80, 0xFF),
		RGB(0xFF, 0xC0, 0xFF),
		RGB(0xC0, 0x00, 0x00),
		RGB(0xC0, 0x40, 0x00),
		RGB(0xC0, 0x80, 0x00),
		RGB(0xC0, 0xC0, 0x00),
		RGB(0xC0, 0xFF, 0x00),
		RGB(0xC0, 0x00, 0x40),
		RGB(0xC0, 0x40, 0x40),
		RGB(0xC0, 0x80, 0x40),
		RGB(0xC0, 0xC0, 0x40),
		RGB(0xC0, 0xFF, 0x40),
		RGB(0xC0, 0x00, 0x80),
		RGB(0xC0, 0x40, 0x80),
		RGB(0xC0, 0x80, 0x80),
		RGB(0xC0, 0xC0, 0x80),
		RGB(0xC0, 0xFF, 0x80),
		RGB(0xC0, 0x00, 0xC0),
		RGB(0xC0, 0x40, 0xC0),
		RGB(0xC0, 0x80, 0xC0),
		RGB(0xC0, 0xFF, 0xC0),
		RGB(0xC0, 0x00, 0xFF),
		RGB(0xC0, 0x40, 0xFF),
		RGB(0xC0, 0x80, 0xFF),
		RGB(0xC0, 0xC0, 0xFF),
		RGB(0xC0, 0xFF, 0xFF),
		RGB(0x80, 0x00, 0x00),
		RGB(0x80, 0x40, 0x00),
		RGB(0x80, 0x80, 0x00),
		RGB(0x80, 0xC0, 0x00),
		RGB(0x80, 0xFF, 0x00),
		RGB(0x80, 0x00, 0x40),
		RGB(0x80, 0x40, 0x40),
		RGB(0x80, 0x80, 0x40),
		RGB(0x80, 0xC0, 0x40),
		RGB(0x80, 0xFF, 0x40),
		RGB(0x80, 0x00, 0x80),
		RGB(0x80, 0x40, 0x80),
		RGB(0x80, 0xC0, 0x80),
		RGB(0x80, 0xFF, 0x80),
		RGB(0x80, 0x00, 0xC0),
		RGB(0x80, 0x40, 0xC0),
		RGB(0x80, 0x80, 0xC0),
		RGB(0x80, 0xC0, 0xC0),
		RGB(0x80, 0xFF, 0xC0),
		RGB(0x80, 0x00, 0xFF),
		RGB(0x80, 0x40, 0xFF),
		RGB(0x80, 0x80, 0xFF),
		RGB(0x80, 0xC0, 0xFF),
		RGB(0x80, 0xFF, 0xFF),
		RGB(0x40, 0x80, 0x00),
		RGB(0x40, 0xC0, 0x00),
		RGB(0x40, 0xFF, 0x00),
		RGB(0x40, 0x80, 0x40),
		RGB(0x40, 0xC0, 0x40),
		RGB(0x40, 0xFF, 0x40),
		RGB(0x40, 0x00, 0x80),
		RGB(0x40, 0x40, 0x80),
		RGB(0x40, 0x80, 0x80),
		RGB(0x40, 0xC0, 0x80),
		RGB(0x40, 0xFF, 0x80),
		RGB(0x40, 0x00, 0xC0),
		RGB(0x40, 0x40, 0xC0),
		RGB(0x40, 0x80, 0xC0),
		RGB(0x40, 0xC0, 0xC0),
		RGB(0x40, 0xFF, 0xC0),
		RGB(0x40, 0x00, 0xFF),
		RGB(0x40, 0x40, 0xFF),
		RGB(0x40, 0x80, 0xFF),
		RGB(0x40, 0xC0, 0xFF),
		RGB(0x40, 0xFF, 0xFF),
		RGB(0x00, 0x80, 0x00),
		RGB(0x00, 0xC0, 0x00),
		RGB(0x00, 0xFF, 0x00),
		RGB(0x00, 0x80, 0x40),
		RGB(0x00, 0xC0, 0x40),
		RGB(0x00, 0xFF, 0x40),
		RGB(0x00, 0x00, 0x80),
		RGB(0x00, 0x40, 0x80),
		RGB(0x00, 0x80, 0x80),
		RGB(0x00, 0xC0, 0x80),
		RGB(0x00, 0xFF, 0x80),
		RGB(0x00, 0x00, 0xC0),
		RGB(0x00, 0x40, 0xC0),
		RGB(0x00, 0x80, 0xC0),
		RGB(0x00, 0xC0, 0xC0),
		RGB(0x00, 0xFF, 0xC0),
		RGB(0x00, 0x00, 0xFF),
		RGB(0x00, 0x40, 0xFF),
		RGB(0x00, 0x80, 0xFF),
		RGB(0x00, 0xC0, 0xFF),
		RGB(0x00, 0xFF, 0xFF)
	};

	m_iColorFPO = (m_iColorFPO + 29) % _countof(colors);
	return colors[m_iColorFPO];
#endif
}

void CFunctionPresentationObject::InitColorIndexes()
{
	m_iColorFPO = 0;
}
