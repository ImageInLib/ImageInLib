
// CursorTestDlg.h : header file
//

#pragma once

#include <vector>

class CCursorButton;

// CCursorTestDlg dialog
class CCursorTestDlg : public CDialogEx
{
// Construction
public:
	CCursorTestDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	enum { IDD = IDD_CURSORTEST_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	int m_nBitDepth = 32;
	std::vector<CCursorButton*> m_cursors;

	void SetBitDepth(int nDepth);

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg void OnDestroy();
	template<int nDepth> afx_msg void OnClickedBitDepth() { SetBitDepth(nDepth); }
	DECLARE_MESSAGE_MAP()
};
