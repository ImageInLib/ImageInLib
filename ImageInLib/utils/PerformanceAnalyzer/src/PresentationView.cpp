
// PresentationView.cpp : implementation of the CPresentationView class
//

#include "stdafx.h"
#include "PerformanceAnalyzer.h"
#include "PresentationView.h"
#include "MemoryDC.h"
#include "Global.h"
#include "ruler.h"
#include "MainFrm.h"
#include "FormatTime.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define ZOOM_STEP 1.7
#define MIN_ZOOM_INTERVAL_WIDTH 3
#define IDTOOLTIP 1

#define CURSOR_COLOR  RGB(255,0,0)
#define MEASURE_MARKER_COLOR  RGB(0,0,255)


namespace
{
	CPoint MakePtInRect(CPoint pt, CRect rc)
	{
		if (pt.x < rc.left)
			pt.x = rc.left;
		if (pt.x > rc.right)
			pt.x = rc.right;
		if (pt.y < rc.top)
			pt.y = rc.top;
		if (pt.y > rc.bottom)
			pt.y = rc.bottom;
		return pt;

	}
}

// CPresentationView

CPresentationView::CPresentationView() :
	m_cursor(Global_C::Instance()->ZoomInfo(), CURSOR_COLOR),
	m_markerMeasureStart(Global_C::Instance()->ZoomInfo(), MEASURE_MARKER_COLOR),
	m_markerMeasureEnd(Global_C::Instance()->ZoomInfo(), MEASURE_MARKER_COLOR),
	m_hBitmap(NULL),
	m_bRedrawCachedObjects(true)
{
	m_mouseAction = NoneAction;

	m_bZooming = false;
	m_bMeasuring = false;
	m_threadListAction = false;
}

CPresentationView::~CPresentationView()
{
}


BEGIN_MESSAGE_MAP(CPresentationView, CWnd)
	ON_WM_PAINT()
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_HSCROLL()
	ON_WM_VSCROLL()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MBUTTONDOWN()
	ON_WM_MBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_CAPTURECHANGED()
	ON_NOTIFY (UDM_TOOLTIP_DISPLAY, NULL, NotifyDisplayTooltip)
	ON_WM_SETCURSOR()
END_MESSAGE_MAP()



// CPresentationView message handlers

BOOL CPresentationView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.style |= WS_VSCROLL | WS_HSCROLL;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	return TRUE;
}

void CPresentationView::NotifyDisplayTooltip(NMHDR * pNMHDR, LRESULT * result)
{
	*result = 0;
	NM_PPTOOLTIP_DISPLAY * pNotify = (NM_PPTOOLTIP_DISPLAY*)pNMHDR;
	if (pNotify->ti->nIDTool != IDTOOLTIP)
		return;
	CPoint pt = *(pNotify->pt);
	ScreenToClient(&pt);
	m_threadList.GetTooltipText(pt, pNotify->ti->sTooltip);	// if return value is 'false' pt is not within threadList object
	                                                        // -> so you can check other objects for tool tip text (ruler,...)
}

int CPresentationView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	int retValue = __super::OnCreate(lpCreateStruct);

	PPTOOLTIP_INFO ti;
	m_tooltip.Create(this, false);
	m_tooltip.SetNotify();
	ti.nBehaviour = PPTOOLTIP_MULTIPLE_SHOW | PPTOOLTIP_DISABLE_AUTOPOP;
	ti.nIDTool = IDTOOLTIP;
	ti.rectBounds = CRect(0, 0, 0, 0);
	ti.sTooltip = "";
	ti.nMask = PPTOOLTIP_MASK_BEHAVIOUR;
	m_tooltip.AddTool(this, ti);

	return retValue;
}

BOOL CPresentationView::PreTranslateMessage(MSG* pMsg)
{
    m_tooltip.RelayEvent(pMsg);
	if (pMsg->message == WM_KEYDOWN && pMsg->wParam == VK_ESCAPE)
		EndMouseAction();
	return __super::PreTranslateMessage(pMsg);
}

bool CPresentationView::CreateCacheBitmap(CDC& dc)
{
	m_bRedrawCachedObjects = true;

	CRect rect;
	GetClientRect(rect);

	if (m_hBitmap)
		return true;
	BITMAPINFO  bmi;
	bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	bmi.bmiHeader.biWidth = rect.Width();
	bmi.bmiHeader.biHeight = rect.Height();
	bmi.bmiHeader.biPlanes = 1;
	bmi.bmiHeader.biBitCount = 24;
	bmi.bmiHeader.biCompression = BI_RGB;
	bmi.bmiHeader.biSizeImage = 0;
	bmi.bmiHeader.biXPelsPerMeter = 0;
	bmi.bmiHeader.biYPelsPerMeter = 0;
	bmi.bmiHeader.biClrUsed = 0;
	bmi.bmiHeader.biClrImportant = 0;
	m_hBitmap = ::CreateDIBSection(dc.m_hDC, &bmi, DIB_RGB_COLORS, NULL, NULL, NULL);
	return (m_hBitmap != nullptr);
}

bool CPresentationView::CreateTextFont(CDC& dc)
{
	LOGFONT lf;
	memset(&lf, 0, sizeof(LOGFONT));
	_tcsncpy_s(lf.lfFaceName, LF_FACESIZE, _T("Segoe UI"), 9);
	lf.lfHeight = -MulDiv(9, GetDeviceCaps(dc.GetSafeHdc(), LOGPIXELSY), 72);
	return (m_font.CreateFontIndirect(&lf) == TRUE);
}

void CPresentationView::DrawCachedObjects(CDC& dc)
{
	CDC* pDC = &dc;

	CRect rectClient;
	GetClientRect(rectClient);

	HBITMAP hOldBmp = nullptr;
	CDC memDC;
	if (memDC.CreateCompatibleDC(pDC))
	{
		hOldBmp = (HBITMAP)memDC.SelectObject(m_hBitmap);
		pDC = &memDC;
	}
	// now draw to bitmap
	if (m_bRedrawCachedObjects)
	{
		// prepare clipping regions
		CRgn rgnRuler, rgnThreadList, rgnClip, rgnClient;
		rgnClip.CreateRectRgn(0, 0, 0, 0);
		rgnRuler.CreateRectRgn(0, 0, rectClient.Width(), 40);	// ruler
		rgnThreadList.CreateRectRgnIndirect(m_threadList.GetDrawingRectangle());
		rgnClip.CombineRgn(&rgnRuler, &rgnThreadList, RGN_OR);
		rgnClient.CreateRectRgnIndirect(rectClient);
		rgnClip.CombineRgn(&rgnClip, &rgnClient, RGN_XOR);

		// draw background
		pDC->SelectClipRgn(&rgnClip);
		CPen pen(PS_SOLID, 1, RGB(255, 255, 255));
		CBrush brush(RGB(255, 255, 255));
		CPen *pOldPen = pDC->SelectObject(&pen);
		CBrush *pOldBrush = pDC->SelectObject(&brush);
		pDC->Rectangle(rectClient);
		pDC->SelectObject(pOldBrush);
		pDC->SelectObject(pOldPen);

		// set text font, color and mode
		CFont* pOldFont = pDC->SelectObject(&m_font);
		pDC->SetTextColor(RGB(0, 0, 0));
		pDC->SetBkMode(TRANSPARENT);

		// draw threads list
		pDC->SelectClipRgn(&rgnThreadList);
		m_threadList.Draw(*pDC);

		// set clip rectangle for ruler
		pDC->SelectClipRgn(&rgnRuler);
		// draw ruler
		DataModel::CZoomInfo & zoomInfo = Global_C::Instance()->ZoomInfo();
		CRect rectRuler(rectClient);
		rectRuler.bottom = rectRuler.top + 40;
		int headerWidth = m_threadList.HeaderWidth();
		rectRuler.left += headerWidth;
		Ruler ruler(*pDC, zoomInfo.ZoomRangeStart(), zoomInfo.ZoomRange(), rectRuler, headerWidth);
		ruler.Draw();
		
		pDC->SelectObject(pOldFont);
	}
	if (memDC.GetSafeHdc() && hOldBmp != nullptr)	// Can copy bitmap?
	{
		dc.BitBlt(0, 0, rectClient.Width(), rectClient.Height(), &memDC, 0, 0, SRCCOPY);
		memDC.SelectObject(hOldBmp);
		if (m_bRedrawCachedObjects)
			m_bRedrawCachedObjects = false;
	}
}

void CPresentationView::DrawCursorObjects(CDC& dc)
{
	CRect cursorRect = m_threadList.GetDataRect();
	cursorRect.top = 0;

	TimeMark markerMeasureStart(m_markerMeasureStart);
	TimeMark markerMeasureEnd(m_markerMeasureEnd);
	if (m_bMeasuring)
	{
		double dStartTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureStart, 0));
		double dEndTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureEnd, 0));
		if (dStartTime > dEndTime)
			std::swap(dStartTime, dEndTime);
		markerMeasureStart.SetTime(dStartTime);
		markerMeasureEnd.SetTime(dEndTime);
	}
	markerMeasureStart.Draw(&dc, cursorRect, TimeMark::LeftArrow);
	markerMeasureEnd.Draw(&dc, cursorRect, TimeMark::RightArrow);
	if (!m_bMeasuring)
		m_cursor.Draw(&dc, cursorRect);

	// set clip rectangle for threads list
	CRgn rgnThreadList;
	rgnThreadList.CreateRectRgnIndirect(m_threadList.GetDrawingRectangle());
	dc.SelectClipRgn(&rgnThreadList);
	if (m_mouseAction == LBAction && m_bZooming)
	{
		CDC tempDC;
		CBitmap bmp;
		CRect rcZoom(m_nZoomStart, m_threadList.GetDataRect().top, m_nZoomEnd, m_threadList.GetDataRect().bottom);
		rcZoom.NormalizeRect();

		tempDC.CreateCompatibleDC(&dc);
		bmp.CreateCompatibleBitmap(&dc, rcZoom.Width(), rcZoom.Height());

		CBitmap* pOldBmp = tempDC.SelectObject(&bmp);
		tempDC.FillSolidRect(0, 0, rcZoom.Width(), rcZoom.Height(), RGB(0, 255, 0));

		BLENDFUNCTION m_bf = { 0 };
		m_bf.BlendOp = AC_SRC_OVER;
		m_bf.SourceConstantAlpha = 64;
		dc.AlphaBlend(rcZoom.left, rcZoom.top,
			rcZoom.Width(), rcZoom.Height(),
			&tempDC,
			0, 0,
			rcZoom.Width(), rcZoom.Height(),
			m_bf);

		tempDC.SelectObject(pOldBmp);
		bmp.DeleteObject();
	}

	if (m_bMeasuring)
	{
		CDC tempDC;
		CBitmap bmp;
		CRect rcMeasure(m_nMeasureStart, m_threadList.GetDataRect().top, m_nMeasureEnd, m_threadList.GetDataRect().bottom);
		rcMeasure.NormalizeRect();

		tempDC.CreateCompatibleDC(&dc);
		bmp.CreateCompatibleBitmap(&dc, rcMeasure.Width(), rcMeasure.Height());

		CBitmap* pOldBmp = tempDC.SelectObject(&bmp);
		tempDC.FillSolidRect(0, 0, rcMeasure.Width(), rcMeasure.Height(), RGB(0, 0, 255));

		BLENDFUNCTION m_bf = { 0 };
		m_bf.BlendOp = AC_SRC_OVER;
		m_bf.SourceConstantAlpha = 64;
		dc.AlphaBlend(rcMeasure.left, rcMeasure.top,
			rcMeasure.Width(), rcMeasure.Height(),
			&tempDC,
			0, 0,
			rcMeasure.Width(), rcMeasure.Height(),
			m_bf);

		tempDC.SelectObject(pOldBmp);
		bmp.DeleteObject();
}
}

void CPresentationView::OnPaint() 
{
	CPaintDC paintDC(this); // device context for painting

	// create font for text drawing
	if (m_font.GetSafeHandle() == nullptr)
		CreateTextFont(paintDC);

	// create bitmap for cache drawing if it does not exist
	if (m_hBitmap == nullptr)
		CreateCacheBitmap(paintDC);

	// create memory dc
	CMemoryDC dc(&paintDC);

	// draw objects
	DrawCachedObjects(dc);
	DrawCursorObjects(dc);
}

void CPresentationView::OnSize(UINT nType, int cx, int cy)
{
	CWnd::OnSize(nType, cx, cy);

	if (m_hBitmap) {
		::DeleteObject(m_hBitmap);
		m_hBitmap = nullptr;
	}

	CRect rectClient;
	GetClientRect(rectClient);
	// set thread list rectangle
	CRect rect(rectClient);
	rectClient.top += 40;	// ruler height
	m_threadList.SetRectangle(rectClient);
	// update scroll bars info
	UpdateScrollInfo();
}

void CPresentationView::UpdateScrollInfo()
{
	SCROLLINFO si;
	m_threadList.GetVScrollInfo(si);
	SetScrollInfo(SB_VERT, &si, TRUE);
	m_threadList.SetVPosition(si.nPos);

	double dPos = m_threadList.GetHScrollInfo(si);
	SetScrollInfo(SB_HORZ, &si, TRUE);
	m_threadList.SetHPosition(dPos);
}

void CPresentationView::UpdateGUI()
{
	UpdateScrollInfo();
	UpdateCursorPosition();
	UpdateZoomRangeTextForStatusBar();
	UpdateMeasuredTimeForStatusBar();
	m_bRedrawCachedObjects = true;
	Invalidate(FALSE);
}

void CPresentationView::ClearView()
{
	m_threadList.Clear();
	m_markerMeasureStart.Clear();
	m_markerMeasureEnd.Clear();
	UpdateGUI();
}

void CPresentationView::InitializeData()
{
	if (m_threadList.Initialize(TPM_Separately))
		UpdateGUI();
}

void CPresentationView::SetPresentationMode(ThreadsPresentationMode presentationMode)
{
	if (PresentationMode() == presentationMode)
		return;
	m_threadList.Clear();
	m_threadList.Initialize(presentationMode);

	UpdateGUI();
}

BOOL CPresentationView::OnEraseBkgnd(CDC*)
{
	return TRUE;
}


void CPresentationView::OnDestroy()
{
	CWnd::OnDestroy();
	if (m_font.GetSafeHandle())
		m_font.DeleteObject();
	if (m_hBitmap) {
		::DeleteObject(m_hBitmap);
		m_hBitmap = nullptr;
	}
}


void CPresentationView::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	// get scroll info
	SCROLLINFO si;
	si.cbSize = sizeof(SCROLLINFO);
	GetScrollInfo(SB_HORZ, &si, SIF_ALL);

	double dMin = 0.;
	double dPage = m_threadList.ScrollInfo_Page();
	double dMax = m_threadList.ScrollInfo_Max() - dPage;
	double dCurPos = m_threadList.ScrollInfo_Pos();
	double dLine = dPage * .05;

	// Determine the new position of scroll box. 
	switch (nSBCode)
	{
		case SB_LEFT:      // Scroll to far left.
			dCurPos = dMin;
			break;

		case SB_RIGHT:      // Scroll to far right.
			dCurPos = dMax;
			break;

		case SB_ENDSCROLL:   // End scroll. 
			break;

		case SB_LINELEFT:      // Scroll left. 
			if (dCurPos > dMin)
				dCurPos = max(dMin, dCurPos - dLine);
			break;

		case SB_LINERIGHT:   // Scroll right. 
			if (dCurPos < dMax)
				dCurPos = min(dMax, dCurPos + dLine);
			break;

		case SB_PAGELEFT:    // Scroll one page left.
		{
			if (dCurPos > dMin)
				dCurPos = max(dMin, dCurPos - dPage);
		}
		break;

		case SB_PAGERIGHT:      // Scroll one page right.
		{
			if (dCurPos < dMax)
				dCurPos = min(dMax, dCurPos + dPage);
		}
		break;

		case SB_THUMBPOSITION:        // Scroll to absolute position. nPos is the position of the scroll box at the end of the drag operation.
		case SB_THUMBTRACK:           // Drag scroll box to specified position. nPos is the position that the scroll box has been dragged to. 
			dCurPos = m_threadList.ScrollInfo_PixelToTime(nPos);
			if (dCurPos < 0.)
				dCurPos = 0.;
			else if (dCurPos > dMax)
				dCurPos = dMax;
			break;
	}

	if (dCurPos != m_threadList.ScrollInfo_Pos()) {
		// set position to scroll bar
		SetScrollPos(SB_HORZ, m_threadList.ScrollInfo_TimeToPixel(dCurPos));
		// set position to thread list
		m_threadList.SetHPosition(dCurPos);
		// update zoom info
		Global_C::Instance()->ZoomInfo().ShiftZoomRangeToPosition(dCurPos);
		UpdateZoomRangeTextForStatusBar();
		// invalidate view
		m_bRedrawCachedObjects = true;
		Invalidate(FALSE);
	}

	CWnd::OnHScroll(nSBCode, nPos, pScrollBar);
}


void CPresentationView::OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar)
{
	// get scroll info
	SCROLLINFO si;
	si.cbSize = sizeof(SCROLLINFO);
	GetScrollInfo(SB_VERT, &si, SIF_ALL);

	int nMinPos = si.nMin;
	int nMaxPos = GetScrollLimit(SB_VERT);
	int nCurPos = si.nPos;
	int nLine = CFunctionPresentationObject::LevelHeight() + CFunctionPresentationObject::LevelSpacing();

	// Determine the new position of scroll box. 
	switch (nSBCode)
	{
		case SB_TOP:      // Scroll to far top.
			nCurPos = nMinPos;
			break;

		case SB_BOTTOM:      // Scroll to far bottom.
			nCurPos = nMaxPos;
			break;

		case SB_ENDSCROLL:   // End scroll. 
			break;

		case SB_LINEUP:      // Scroll up. 
			if (nCurPos > nMinPos)
				nCurPos = max(nMinPos, nCurPos - nLine);
			break;

		case SB_LINEDOWN:   // Scroll down. 
			if (nCurPos < nMaxPos)
				nCurPos = min(nMaxPos, nCurPos + nLine);
			break;

		case SB_PAGEUP:    // Scroll one page up.
		{
			if (nCurPos > nMinPos)
				nCurPos = max(nMinPos, nCurPos - (int)si.nPage);
		}
		break;

		case SB_PAGEDOWN:      // Scroll one page down.
		{
			if (nCurPos < nMaxPos)
				nCurPos = min(nMaxPos, nCurPos + (int)si.nPage);
		}
		break;

		case SB_THUMBPOSITION:        // Scroll to absolute position. nPos is the position
			nCurPos = si.nTrackPos;     // of the scroll box at the end of the drag operation. 
			break;

		case SB_THUMBTRACK:           // Drag scroll box to specified position. nPos is the
			nCurPos = si.nTrackPos;     // position that the scroll box has been dragged to. 
			break;
	}

	if (si.nPos != nCurPos) {
		SetScrollPos(SB_VERT, nCurPos);
		m_threadList.SetVPosition(nCurPos);
		// invalidate view
		m_bRedrawCachedObjects = true;
		Invalidate(FALSE);
	}

	CWnd::OnVScroll(nSBCode, nPos, pScrollBar);
}

bool CPresentationView::StartMouseAction(MouseButtonAction action)
{
	if (m_mouseAction != NoneAction)
		return false;

	m_mouseAction = action;
	SetCapture();
	SetFocus();
	return true;
}

void CPresentationView::EndMouseAction()
{
	if (m_mouseAction == NoneAction)
		return;
	
	m_mouseAction = NoneAction;
	m_bZooming = false;
	m_bMeasuring = false;
	ReleaseCapture();
}

void CPresentationView::OnLButtonDown(UINT nFlags, CPoint point)
{
	if (StartMouseAction(LBAction))
	{
		CRect rect = m_threadList.GetDataRect();
		if (rect.PtInRect(point))
		{
			m_threadListAction = true;
			m_nZoomStart = m_nZoomEnd = point.x;
		}
		else
		{
			StartMeasurement(point);
		}
	}
	CWnd::OnLButtonDown(nFlags, point);
}

void CPresentationView::OnLButtonUp(UINT nFlags, CPoint point)
{
	if (m_mouseAction != LBAction)
		return CWnd::OnLButtonUp(nFlags, point);

	if (m_bZooming)
	{
		m_nZoomEnd = MakePtInRect(point, m_threadList.GetDataRect()).x;
		
		if (abs(m_nZoomEnd - m_nZoomStart) <= MIN_ZOOM_INTERVAL_WIDTH )
			m_bZooming = false;
	}

	if (m_bZooming)
	{
		// zooming
		double dStartTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nZoomStart, 0));
		double dEndTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nZoomEnd, 0));
		if (dStartTime > dEndTime)
			std::swap(dStartTime, dEndTime);

		ZoomToTimeInterval(dStartTime, dEndTime);
		m_bZooming = false;
		Invalidate(0);
	}
	else if (m_bMeasuring)
	{
		EndMeasurement(point);
	}
	else
	{
		if (m_threadListAction && m_threadList.GetDataRect().PtInRect(point))
			ZoomInOut(m_threadList.GetRelativeTimeFromPixel(point), ZOOM_STEP);
	}

	EndMouseAction();
	CWnd::OnLButtonUp(nFlags, point);
}

void CPresentationView::OnMButtonDown(UINT nFlags, CPoint point)
{
	if (!StartMouseAction(MBAction))
		return CWnd::OnLButtonDown(nFlags, point);

	if (m_threadList.GetDataRect().PtInRect(point))
	{
		m_threadListAction = true;
	}

	CWnd::OnMButtonDown(nFlags, point);
}

void CPresentationView::OnMButtonUp(UINT nFlags, CPoint point)
{
	if (m_mouseAction != MBAction)
		return CWnd::OnLButtonUp(nFlags, point);

	if (m_threadListAction && m_threadList.GetDataRect().PtInRect(point))
	{
		ZoomToTimeInterval(0.0, Global_C::Instance()->ZoomInfo().DataWidth());
	}

	EndMouseAction();
	CWnd::OnMButtonUp(nFlags, point);
}

void CPresentationView::OnRButtonDown(UINT nFlags, CPoint point)
{
	if (StartMouseAction(RBAction))
		StartMeasurement(point);
	CWnd::OnRButtonDown(nFlags, point);
}

void CPresentationView::OnRButtonUp(UINT nFlags, CPoint point)
{
	if (m_mouseAction == RBAction)
	{
		if (m_bMeasuring)
		{
			EndMeasurement(point);
		}
		else
		{
			if (m_threadListAction && m_threadList.GetDataRect().PtInRect(point))
				ZoomInOut(m_threadList.GetRelativeTimeFromPixel(point), 1/ZOOM_STEP);
		}

		EndMouseAction();
	}
	CWnd::OnRButtonUp(nFlags, point);
}

void CPresentationView::OnMouseMove(UINT nFlags, CPoint point)
{
	UpdateCursorPosition();

	if (m_threadListAction && m_mouseAction == LBAction && !m_bMeasuring)
	{
		m_nZoomEnd = MakePtInRect(point, m_threadList.GetDataRect()).x;
		if (abs(m_nZoomEnd - m_nZoomStart) > MIN_ZOOM_INTERVAL_WIDTH)
			m_bZooming = true;
		else
			m_bZooming = false;
		Invalidate(0);
	}

	if (m_threadListAction && (m_mouseAction == RBAction || m_bMeasuring))
	{
		int pos = MakePtInRect(point, m_threadList.GetDataRect()).x;
		if (!m_bMeasuring)
		{
			if (pos != m_nMeasureStart)
			{
				m_bMeasuring = true;
				m_nMeasureEnd = pos;
			}
		}
		else
		{
			if (m_bMeasureMovingEnd)
				m_nMeasureEnd = pos;
			else
				m_nMeasureStart = pos;
		}
		Invalidate(0);
		UpdateMeasuredTimeForStatusBar();
	}

	CWnd::OnMouseMove(nFlags, point);
}

void CPresentationView::OnCaptureChanged(CWnd *pWnd)
{
	m_mouseAction = NoneAction;

	m_threadListAction = false;
	m_bZooming = false;
	m_bMeasuring = false;
	Invalidate(0);

	CWnd::OnCaptureChanged(pWnd);
}

void CPresentationView::UpdateCursorPosition()
{
	if (Global_C::Instance()->ZoomInfo().DataWidth() <= 0.0)
		return;

	CString strTimePos;
	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(&pt);	
	pt = MakePtInRect(pt, m_threadList.GetDataRect());
	double dTime = m_threadList.GetRelativeTimeFromPixel(pt);
	if (m_cursor.SetTime(dTime))
		Invalidate(0);

	strTimePos = FormatAbsoluteTime(m_cursor.GetTime() + Global_C::Instance()->ZoomInfo().DataShift(), DefaultPrecisionAfterSeconds);
	((CMainFrame*)theApp.GetMainWnd())->StatusBarText(ID_STATUSBAR_TIMEPOSITION) = strTimePos;

}

void CPresentationView::ZoomInOut(double dCenterTime, double dZoomFactor)
{
	DataModel::CZoomInfo& zoomInfo = Global_C::Instance()->ZoomInfo();
	double dZoomRange = min(zoomInfo.ZoomRange() / dZoomFactor, zoomInfo.DataWidth());

	double dStart = dCenterTime - 0.5*dZoomRange;
	if (dStart < 0.0)
		dStart = 0.0;
	if (dStart > zoomInfo.DataWidth() - dZoomRange)
		dStart = zoomInfo.DataWidth() - dZoomRange;

	ZoomToTimeInterval(dStart, dStart + dZoomRange);

}

void CPresentationView::ZoomToTimeInterval(double dStart, double dEnd)
{
	DataModel::CZoomInfo& zoomInfo = Global_C::Instance()->ZoomInfo();

	double dRange = max(dEnd - dStart, zoomInfo.ZoomMinRange());
	if (dRange == zoomInfo.ZoomRange())
		return;

	zoomInfo.SetZoomRange(dStart, dRange);
	m_threadList.ZoomInfoChanged();
	m_threadList.SetHPosition(dStart);
	UpdateGUI();
}

void CPresentationView::UpdateZoomRangeTextForStatusBar()
{
	DataModel::CZoomInfo& zoomInfo = Global_C::Instance()->ZoomInfo();

	((CMainFrame*)theApp.GetMainWnd())->StatusBarText(ID_STATUSBAR_TIMERANGE) =
		FormatTimeInterval(zoomInfo.ZoomRange(), DefaultIntervalPrecision);

}

void CPresentationView::UpdateMeasuredTimeForStatusBar()
{
	TimeMark markerMeasureStart(m_markerMeasureStart);
	TimeMark markerMeasureEnd(m_markerMeasureEnd);
	if (m_bMeasuring)
	{
		double dStartTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureStart, 0));
		double dEndTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureEnd, 0));
		if (dStartTime > dEndTime)
			std::swap(dStartTime, dEndTime);
		markerMeasureStart.SetTime(dStartTime);
		markerMeasureEnd.SetTime(dEndTime);
	}
	if (markerMeasureStart.IsValid() && markerMeasureEnd.IsValid())
		((CMainFrame*)theApp.GetMainWnd())->StatusBarText(ID_STATUSBAR_MEASUREDTIME) =
			FormatTimeInterval(markerMeasureEnd.GetTime() - markerMeasureStart.GetTime(), DefaultIntervalPrecision);
	else
		((CMainFrame*)theApp.GetMainWnd())->StatusBarText(ID_STATUSBAR_MEASUREDTIME) =	L"";
}

void CPresentationView::SetMeasureMarkers(double dStart, double dEnd)
{
	bool bUpdated = false;
	bUpdated |= m_markerMeasureStart.SetTime(dStart);
	bUpdated |= m_markerMeasureEnd.SetTime(dEnd);
	if (bUpdated)
	{
		Invalidate(0);
		UpdateMeasuredTimeForStatusBar();
	}
}

int CPresentationView::HitTestMeasureMarkers(CPoint point)
{
	CRect rect = m_threadList.GetDataRect();
	rect.bottom = rect.top;
	rect.top = 0;
	if (rect.PtInRect(point))
	{
		if (m_markerMeasureStart.IsValid() && m_markerMeasureEnd.IsValid())
		{
			int start = m_threadList.GetPixelPositionFromRelativeTime(m_markerMeasureStart.GetTime());
			if (abs(point.x - start) <= 2)
				return 0;
		}
		if (m_markerMeasureEnd.IsValid())
		{
			int end = m_threadList.GetPixelPositionFromRelativeTime(m_markerMeasureEnd.GetTime());
			if (abs(point.x - end) <= 2)
				return 1;
		}
	}
	return -1;
}

bool CPresentationView::StartMeasurement(CPoint point)
{
	CRect rect = m_threadList.GetDataRect();
	if (rect.PtInRect(point))
	{
		m_threadListAction = true;
		m_nMeasureStart = m_nMeasureEnd = point.x;
		m_bMeasureMovingEnd = true;
		return true;
	}

	if (m_markerMeasureStart.IsValid() && m_markerMeasureEnd.IsValid())
	{
		int hit = HitTestMeasureMarkers(point);
		if (hit == 0 || hit == 1)
		{
			m_bMeasureMovingEnd = (hit == 1);
			m_nMeasureStart = m_threadList.GetPixelPositionFromRelativeTime(m_markerMeasureStart.GetTime());
			m_nMeasureEnd = m_threadList.GetPixelPositionFromRelativeTime(m_markerMeasureEnd.GetTime());
			m_threadListAction = true;
			m_bMeasuring = true;
			Invalidate(0);
			return true;
		}
	}

	return false;
}

void CPresentationView::EndMeasurement(CPoint point)
{
	int pos = MakePtInRect(point, m_threadList.GetDataRect()).x;
	if (m_bMeasureMovingEnd)
		m_nMeasureEnd = pos;
	else
		m_nMeasureStart = pos;

	// measuring time
	double dStartTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureStart, 0));
	double dEndTime = m_threadList.GetRelativeTimeFromPixel(CPoint(m_nMeasureEnd, 0));
	if (dStartTime > dEndTime)
		std::swap(dStartTime, dEndTime);

	SetMeasureMarkers(dStartTime, dEndTime);
	m_bMeasuring = false;
	Invalidate(0);
}

BOOL CPresentationView::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message)
{
	POINT pt;
	GetCursorPos(&pt);
	ScreenToClient(&pt);	
	int hit = HitTestMeasureMarkers(pt);
	if (hit == 0 || hit == 1)
	{
		::SetCursor(::LoadCursor(NULL, IDC_SIZEWE));
		return TRUE;
	}
	return __super::OnSetCursor(pWnd, nHitTest, message);
}
