
// CursorTestDlg.cpp : implementation file
//

#include "stdafx.h"
#include "CursorTest.h"
#include "CursorTestDlg.h"
#include "afxdialogex.h"
#include "ReducedColorsCursor.h"

class CCursorButton : public CButton
{
public:
	CCursorButton(UINT nCursorID) : CButton(), m_nCursorID(nCursorID) {}
	virtual ~CCursorButton() {}

	void SetBitDepth(int nBitDepth)
	{
		m_nBitDepth = nBitDepth;
		m_hCursor = nullptr;
		Invalidate();
	}

	void DrawItem(LPDRAWITEMSTRUCT lpDrawItemStruct) override
	{
		UINT uStyle = DFCS_BUTTONPUSH;

		// This code only works with buttons.
		ASSERT(lpDrawItemStruct->CtlType == ODT_BUTTON);

		// If drawing selected, add the pushed style to DrawFrameControl.
		if (m_bDown || (lpDrawItemStruct->itemState & ODS_SELECTED))
			uStyle |= DFCS_PUSHED;

		// Draw the button frame.
		::DrawFrameControl(lpDrawItemStruct->hDC, &lpDrawItemStruct->rcItem, DFC_BUTTON, uStyle);

		if (!m_hCursor)
		{
			m_hCursor = ReducedColorsCursor::LoadCursor(m_nCursorID, m_nBitDepth);

			m_size = { GetSystemMetrics(SM_CXCURSOR), GetSystemMetrics(SM_CYCURSOR) };
			ICONINFO info = { 0 };
			if (::GetIconInfo(m_hCursor, &info))
			{
				bool bBWCursor = !info.hbmColor;
				BITMAP bmpinfo = { 0 };
				if (::GetObject(info.hbmMask, sizeof(BITMAP), &bmpinfo))
				{
					m_size = { bmpinfo.bmWidth, abs(bmpinfo.bmHeight) / (bBWCursor ? 2 : 1) };
				}

				::DeleteObject(info.hbmColor);
				::DeleteObject(info.hbmMask);
			}
		}

		int xOffset = (lpDrawItemStruct->rcItem.right - lpDrawItemStruct->rcItem.left - m_size.cx) / 2;
		int yOffset = (lpDrawItemStruct->rcItem.bottom - lpDrawItemStruct->rcItem.top - m_size.cy) / 2;

		::DrawIcon(lpDrawItemStruct->hDC, xOffset, yOffset, m_hCursor);
	}

protected:
	DECLARE_MESSAGE_MAP()

	afx_msg void OnLButtonDown(UINT nFlags, CPoint point)
	{
		GetCapture();
		m_bDown = true;
		::SetCursor(m_hCursor);
		__super::OnLButtonDown(nFlags, point);
	}

	afx_msg void OnLButtonUp(UINT nFlags, CPoint point)
	{
		m_bDown = false;
		ReleaseCapture();
		__super::OnLButtonUp(nFlags, point);
		Invalidate();
	}

	afx_msg void OnCaptureChanged(CWnd * pWnd)
	{
		m_bDown = false;
		__super::OnCaptureChanged(pWnd);
		Invalidate();
	}

protected:
	int m_nBitDepth { 32 };
	UINT m_nCursorID { (UINT)-1 };
	HCURSOR m_hCursor { nullptr };
	CSize m_size { 0,0 };
	bool m_bDown { false };
};

BEGIN_MESSAGE_MAP(CCursorButton, CButton)
	ON_WM_CAPTURECHANGED()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
END_MESSAGE_MAP()

// CCursorTestDlg dialog

CCursorTestDlg::CCursorTestDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CCursorTestDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CCursorTestDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

void CCursorTestDlg::SetBitDepth(int nDepth)
{
	m_nBitDepth = nDepth;

	for (auto ptrCursor : m_cursors)
	{
		ptrCursor->SetBitDepth(nDepth);
	}
}

BEGIN_MESSAGE_MAP(CCursorTestDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_RADIO_1BPP, OnClickedBitDepth<1>)
	ON_BN_CLICKED(IDC_RADIO_4BPP, OnClickedBitDepth<4>)
	ON_BN_CLICKED(IDC_RADIO_8BPP, OnClickedBitDepth<8>)
	ON_BN_CLICKED(IDC_RADIO_24BPP, OnClickedBitDepth<24>)
	ON_BN_CLICKED(IDC_RADIO_32BPP, OnClickedBitDepth<32>)
	ON_WM_DESTROY()
END_MESSAGE_MAP()


// CCursorTestDlg message handlers

void CCursorTestDlg::OnDestroy()
{
	for (auto ptr : m_cursors)
	{
		ptr->DestroyWindow();
		::delete ptr;
	}

	m_cursors.clear();

	__super::OnDestroy();
}

BOOL CCursorTestDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	((CButton*)GetDlgItem(IDC_RADIO_32BPP))->SetCheck(true);


	const int size = 48;
	const int delim = 6;

	CRect rct;
	GetClientRect(rct);

	rct.DeflateRect({ 2*delim, 2*delim });
	rct.left += 80;
	rct.right -= size;

	CPoint ptPos = rct.TopLeft();

	for (auto nID : { IDC_CURSOR, IDC_CROSS, IDC_CROSS_BW })
	{
		auto ptr = ::new CCursorButton(nID);
		VERIFY(ptr);

		ptr->Create(L"", WS_VISIBLE | WS_CHILD | BS_OWNERDRAW, CRect { ptPos, CSize {size,size} }, this, nID);

		m_cursors.push_back(ptr);

		ptPos.x += size + delim;
		if (ptPos.x > rct.right)
		{
			ptPos.y += size + delim;
			ptPos.x = rct.left;
		}
	}

	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CCursorTestDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CCursorTestDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

