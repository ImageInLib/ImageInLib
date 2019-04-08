#pragma once

#include "ThreadsList.h"
#include "TimeMark.h"
#include "PPTooltip\PPTooltip.h"


enum MouseButtonAction
{
	NoneAction,
	LBAction,
	MBAction,
	RBAction
};

class CPresentationView : public CWnd
{
public:
	CPresentationView();
	virtual ~CPresentationView();

	void ClearView();
	void InitializeData();

	ThreadsPresentationMode PresentationMode() const { return m_threadList.PresentationMode(); }
	void SetPresentationMode(ThreadsPresentationMode presentationMode);

	void UpdateZoomRangeTextForStatusBar();
	void UpdateMeasuredTimeForStatusBar();

protected:
	void UpdateScrollInfo();
	void UpdateGUI();	// update cursor position, zoom range text and scroll info

private:
	bool CreateCacheBitmap(CDC& dc);
	bool CreateTextFont(CDC& dc);
	void DrawCachedObjects(CDC& dc);
	void DrawCursorObjects(CDC& dc);

protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual BOOL PreTranslateMessage(MSG* pMsg);
	DECLARE_MESSAGE_MAP()
	afx_msg void OnPaint();
	afx_msg void OnDestroy();
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

public:
	afx_msg void OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);
	afx_msg void OnVScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar);

	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnMButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnCaptureChanged(CWnd *pWnd);
	afx_msg void NotifyDisplayTooltip(NMHDR * pNMHDR, LRESULT * result);
	afx_msg	BOOL OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message);

	void ZoomInOut(double dCenterTime, double dZoomStep);
	void ZoomToTimeInterval(double dStart, double dEnd);
	void UpdateCursorPosition();
	bool StartMouseAction(MouseButtonAction action);
	void EndMouseAction();
	void SetMeasureMarkers(double dStart, double dEnd);
	bool StartMeasurement(CPoint point);
	void EndMeasurement(CPoint point);
	int HitTestMeasureMarkers(CPoint point);

protected:
	CPPToolTip		m_tooltip;
	CThreadsList	m_threadList;
	TimeMark		m_cursor;
	TimeMark		m_markerMeasureStart;
	TimeMark		m_markerMeasureEnd;
	CFont			m_font;
	bool			m_bZooming;
	int				m_nZoomStart;
	int				m_nZoomEnd;
	bool			m_bMeasuring;
	bool			m_bMeasureMovingEnd;
	int				m_nMeasureStart;
	int				m_nMeasureEnd;
	HBITMAP   m_hBitmap;
	bool      m_bRedrawCachedObjects;
	
	MouseButtonAction m_mouseAction;

	bool			m_threadListAction;
};
