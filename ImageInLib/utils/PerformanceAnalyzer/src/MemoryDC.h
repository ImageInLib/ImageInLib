#pragma once

//////////////////////////////////////////////////
// CMemoryDC - memory DC
//
// Author: Keith Rule
// Email:  keithr@europa.com
// Copyright 1996-1997, Keith Rule
//
// You may freely use or modify this code provided this
// Copyright is included in all derived versions.
//
// History - 10/3/97 Fixed scrolling bug.
//                   Added print support.
//           25 feb 98 - fixed minor assertion bug
//
// This class implements a memory Device Context

class CMemoryDC : public CDC
{
	//Attributes
private:
	HBITMAP		m_hOldBmp;

	CDC*		m_pDC;         // Saves CDC passed in constructor
	CRect		m_rect;        // Rectangle of drawing area.
	BOOL		m_bMemDC;      // TRUE if CDC really is a Memory DC.

	HBITMAP		m_hBitmap;
	BITMAPINFO  m_bmi;

	//Public interface methods
public:
	const CRect &DrawingRect(void) const { return (m_rect); }
	void SetDrawingRect(const CRect &r) { m_rect = r; }

	// constructor sets up the memory DC
	CMemoryDC(CDC* pDC, HBITMAP hBitmap = NULL) : CDC()
	{
		ASSERT(pDC != NULL);

		m_pDC = pDC;
		m_bMemDC = !pDC->IsPrinting();

		if (m_bMemDC && !CreateCompatibleDC(pDC)) {
			m_bMemDC = false;
			return;
		}

		m_hBitmap = hBitmap;
		if (m_bMemDC)    // Create a Memory DC
		{
			pDC->GetClipBox(&m_rect);

			HBITMAP hMemBmp;
			if (!m_hBitmap) {
				m_bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
				m_bmi.bmiHeader.biWidth = m_rect.Width();
				m_bmi.bmiHeader.biHeight = m_rect.Height();
				m_bmi.bmiHeader.biPlanes = 1;
				m_bmi.bmiHeader.biBitCount = 24;
				m_bmi.bmiHeader.biCompression = BI_RGB;
				m_bmi.bmiHeader.biSizeImage = 0;
				m_bmi.bmiHeader.biXPelsPerMeter = 0;
				m_bmi.bmiHeader.biYPelsPerMeter = 0;
				m_bmi.bmiHeader.biClrUsed = 0;
				m_bmi.bmiHeader.biClrImportant = 0;
				hMemBmp = ::CreateDIBSection(pDC->m_hDC, &m_bmi, DIB_RGB_COLORS, NULL, NULL, NULL);
				//::CreateDIBitmap(pDC->m_hDC, (BITMAPINFOHEADER *) pBackBmpInfo, CBM_INIT, pBackBmpData, pBackBmpInfo, DIB_RGB_COLORS);
			}
			else {
				hMemBmp = m_hBitmap;
			}
			m_hOldBmp = (HBITMAP) ::SelectObject(m_hDC, hMemBmp);

			//+++HBITMAP hMemBmp = ::CreateCompatibleBitmap(pDC->m_hDC, m_rect.Width(), m_rect.Height());
			//+++m_hOldBmp = (HBITMAP) ::SelectObject(m_hDC, hMemBmp);

			SetWindowOrg(m_rect.left, m_rect.top);
		}
		else        // Make a copy of the relevent parts of the current DC for printing
		{
			m_bPrinting = pDC->m_bPrinting;
			m_hDC = pDC->m_hDC;
			m_hAttribDC = pDC->m_hAttribDC;
		}
	}

	// Destructor copies the contents of the mem DC to the original DC
	~CMemoryDC()
	{
		if (m_bMemDC)
		{
			// Copy the offscreen bitmap onto the screen.
			m_pDC->BitBlt(m_rect.left, m_rect.top, m_rect.Width(), m_rect.Height(), this, m_rect.left, m_rect.top, SRCCOPY);

			//Swap back the original bitmap.
			HBITMAP hBitmap = (HBITMAP)::SelectObject(m_hDC, m_hOldBmp);
			if (!m_hBitmap)
				::DeleteObject(hBitmap);

		}
		else {
			// All we need to do is replace the DC with an illegal value,
			// this keeps us from accidently deleting the handles associated with
			// the CDC that was passed to the constructor.
			m_hDC = m_hAttribDC = NULL;
		}
	}

	// Allow usage as a pointer
	CMemoryDC* operator->() { return this; }

	// Allow usage as a pointer
	operator CMemoryDC*() { return this; }
};
