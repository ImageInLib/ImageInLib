
#pragma once

#include "PresentationView.h"

class CMainFrame : public CFrameWnd
{
public:
	CMainFrame();
	virtual ~CMainFrame();
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	virtual BOOL OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo);
	CString& StatusBarText(UINT nID){ return m_statusBarString[nID]; }

#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:  // control bar embedded members
	CToolBar          m_wndToolBar;
	CStatusBar        m_wndStatusBar;
	CPresentationView m_wndView;

	std::map<UINT, CString> m_statusBarString;

protected:
	DECLARE_DYNAMIC(CMainFrame)
	DECLARE_MESSAGE_MAP()

	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSetFocus(CWnd *pOldWnd);
	afx_msg void OnApplicationLook(UINT id);
	afx_msg void OnUpdateApplicationLook(CCmdUI* pCmdUI);
	afx_msg void OnOpenLog();
	afx_msg void OnUpdateOpenLog(CCmdUI* pCmdUI);
	afx_msg void OnUpdateMessage(CCmdUI *pCmdUI);
	afx_msg void OnUpdateCursorPosition(CCmdUI *pCmdUI);
	afx_msg void OnUpdateMeasuredTime(CCmdUI *pCmdUI);
	afx_msg void OnUpdateTimeRange(CCmdUI *pCmdUI);
	afx_msg void OnViewDisplaythreadsseparately();
	afx_msg void OnUpdateViewDisplaythreadsseparately(CCmdUI *pCmdUI);
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);
	virtual BOOL PreTranslateMessage(MSG * pMsg);

	bool HandleScrollbars(WPARAM wParam);
};


