
// PerformanceAnalyzer.h : main header file for the PerformanceAnalyzer application
//
#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"       // main symbols


// CPerformanceAnalyzerApp:
// See PerformanceAnalyzer.cpp for the implementation of this class
//

class CPerformanceAnalyzerApp : public CWinApp
{
public:
	CPerformanceAnalyzerApp();


// Overrides
public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();

// Implementation

public:
	UINT  m_nAppLook;
	afx_msg void OnAppAbout();
	DECLARE_MESSAGE_MAP()
};

extern CPerformanceAnalyzerApp theApp;
