
// MainFrm.cpp : implementation of the CMainFrame class
//

#include "stdafx.h"
#include "PerformanceAnalyzer.h"
#include "MainFrm.h"
#include "Global.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CMainFrame

IMPLEMENT_DYNAMIC(CMainFrame, CFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
	ON_WM_CREATE()
	ON_WM_SETFOCUS()
	ON_WM_GETMINMAXINFO()
	ON_COMMAND(ID_OPEN_LOG, &CMainFrame::OnOpenLog)
	ON_UPDATE_COMMAND_UI(ID_OPEN_LOG, &CMainFrame::OnUpdateOpenLog)
	ON_COMMAND(ID_VIEW_DISPLAYTHREADSSEPARATELY, &CMainFrame::OnViewDisplaythreadsseparately)
	ON_UPDATE_COMMAND_UI(ID_VIEW_DISPLAYTHREADSSEPARATELY, &CMainFrame::OnUpdateViewDisplaythreadsseparately)
	ON_UPDATE_COMMAND_UI(ID_STATUSBAR_MESSAGE, OnUpdateMessage)
	ON_UPDATE_COMMAND_UI(ID_STATUSBAR_TIMEPOSITION, OnUpdateCursorPosition)
	ON_UPDATE_COMMAND_UI(ID_STATUSBAR_MEASUREDTIME, OnUpdateMeasuredTime)
	ON_UPDATE_COMMAND_UI(ID_STATUSBAR_TIMERANGE, OnUpdateTimeRange)

END_MESSAGE_MAP()

// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
	// TODO: add member initialization code here
}

CMainFrame::~CMainFrame()
{
}

int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CFrameWnd::OnCreate(lpCreateStruct) == -1)
		return -1;

	// create a view to occupy the client area of the frame
	if (!m_wndView.Create(NULL, NULL, AFX_WS_DEFAULT_VIEW, CRect(0, 0, 0, 0), this, AFX_IDW_PANE_FIRST, NULL))
	{
		TRACE0("Failed to create view window\n");
		return -1;
	}

	if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | CBRS_TOP | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC ) ||
		!m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}

	if (!m_wndStatusBar.Create(this))
	{
		TRACE0("Failed to create status bar\n");
		return -1;      // fail to create
	}

	UINT indicators[] = { ID_STATUSBAR_MESSAGE, ID_STATUSBAR_TIMEPOSITION, ID_STATUSBAR_MEASUREDTIME };
	m_wndStatusBar.SetIndicators(indicators, _countof(indicators));
	m_wndStatusBar.SetPaneInfo(0, ID_STATUSBAR_MESSAGE, SBPS_STRETCH, 0);	
	m_wndStatusBar.SetPaneInfo(1, ID_STATUSBAR_TIMEPOSITION, 0, 80);	
	m_wndStatusBar.SetPaneInfo(2, ID_STATUSBAR_MEASUREDTIME, 0, 80);	

	// TODO: Delete these three lines if you don't want the toolbar to be dockable
	m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
	EnableDocking(CBRS_ALIGN_ANY);
	DockControlBar(&m_wndToolBar);

	return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
	if( !CFrameWnd::PreCreateWindow(cs) )
		return FALSE;
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	cs.style = WS_OVERLAPPED | WS_CAPTION | FWS_ADDTOTITLE
		 | WS_THICKFRAME | WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_MAXIMIZE | WS_SYSMENU;

	cs.dwExStyle &= ~WS_EX_CLIENTEDGE;
	cs.lpszClass = AfxRegisterWndClass(0);
	return TRUE;
}

// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
	CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
	CFrameWnd::Dump(dc);
}
#endif //_DEBUG


// CMainFrame message handlers

void CMainFrame::OnSetFocus(CWnd* /*pOldWnd*/)
{
	// forward focus to the view window
	m_wndView.SetFocus();
}

BOOL CMainFrame::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
	// let the view have first crack at the command
	if (m_wndView.OnCmdMsg(nID, nCode, pExtra, pHandlerInfo))
		return TRUE;

	// otherwise, do default handling
	return CFrameWnd::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

void CMainFrame::OnOpenLog()
{
	CFileDialog fileDialog(TRUE, _T("txt"), NULL, OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST, NULL, this, 0, TRUE);
	if (fileDialog.DoModal() != IDOK)
		return;

	// clear presentation view and old data
	CWaitCursor wait;
	m_wndView.ClearView();
	Global_C::Instance()->ThreadContainer().Clear();
	Global_C::Instance()->ThreadContainer().ProcessLogFile(fileDialog.GetPathName());
	m_wndView.InitializeData();
}

void CMainFrame::OnUpdateOpenLog(CCmdUI* pCmdUI)
{
	pCmdUI->Enable(TRUE);
}

void CMainFrame::OnViewDisplaythreadsseparately()
{
	ThreadsPresentationMode presentationMode = m_wndView.PresentationMode() == TPM_Separately ? TPM_Jointly : TPM_Separately;
	m_wndView.SetPresentationMode(presentationMode);
}

void CMainFrame::OnUpdateViewDisplaythreadsseparately(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_wndView.PresentationMode() == TPM_Separately);
}


void CMainFrame::OnGetMinMaxInfo(MINMAXINFO* lpMMI)
{
	//lpMMI->ptMinTrackSize.x = 400;
	//lpMMI->ptMinTrackSize.y = 400;

	CFrameWnd::OnGetMinMaxInfo(lpMMI);
}

void CMainFrame::OnUpdateMeasuredTime(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetText(StatusBarText(ID_STATUSBAR_MEASUREDTIME));
}

void CMainFrame::OnUpdateTimeRange(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetText(StatusBarText(ID_STATUSBAR_TIMERANGE));
}

void CMainFrame::OnUpdateMessage(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetText(StatusBarText(ID_STATUSBAR_MESSAGE));
}

void CMainFrame::OnUpdateCursorPosition(CCmdUI *pCmdUI)
{
	pCmdUI->Enable();
	pCmdUI->SetText(StatusBarText(ID_STATUSBAR_TIMEPOSITION));
}

BOOL CMainFrame::PreTranslateMessage(MSG * pMsg)
{
	if (pMsg->message == WM_KEYDOWN)
	{
		if (HandleScrollbars(pMsg->wParam))
			return TRUE;
	}

	if (pMsg->message == WM_MOUSEWHEEL)
	{
		bool bPositive = ((short)HIWORD(pMsg->wParam) > 0);
		m_wndView.SendMessage(WM_HSCROLL, MAKELONG(bPositive ? SB_LINEUP : SB_LINEDOWN, 0), 0L);
		return TRUE;
	}

	return __super::PreTranslateMessage(pMsg);
}

bool CMainFrame::HandleScrollbars(WPARAM wParam)
{
	USHORT wScrollNotify = 0;
	UINT uMessage = 0;
	switch (wParam) 
    {
        case VK_UP:
            wScrollNotify = SB_LINEUP;
            uMessage = WM_VSCROLL;
            break;

		case VK_DOWN:
            wScrollNotify = SB_LINEDOWN;
            uMessage = WM_VSCROLL;
            break;

        case VK_PRIOR:    //PAGEUP key
            wScrollNotify = SB_PAGEDOWN;
            uMessage = WM_HSCROLL;
            break;

        case VK_NEXT:     // PAGEDOWN key
            wScrollNotify = SB_PAGEUP;
            uMessage = WM_HSCROLL;
            break;

        case VK_HOME:
            wScrollNotify = SB_BOTTOM;
            uMessage = WM_HSCROLL;
            break;

        case VK_END:
            wScrollNotify = SB_TOP;
            uMessage = WM_HSCROLL;
            break;

        case VK_RIGHT:
            wScrollNotify = SB_LINEDOWN;
            uMessage = WM_HSCROLL;
            break;

        case VK_LEFT:
            wScrollNotify = SB_LINEUP;
            uMessage = WM_HSCROLL;
            break;

        default:
			return false;
     }

     m_wndView.SendMessage(uMessage, MAKELONG(wScrollNotify, 0), 0L);
	 return true;
}
