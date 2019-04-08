#pragma once

#include <vector>
#include "Thread.h"
#include "ThreadPresentationObject.h"

enum ThreadsPresentationMode { TPM_Separately, TPM_Jointly };

class CThreadsListFiller;

class CThreadsList
{
	friend CThreadsListFiller;
public:
	CThreadsList();
	~CThreadsList();

	bool Initialize(ThreadsPresentationMode presentationMode);
	void Clear();
	void Draw(CDC &dc);

	LONG HeaderWidth() const;
	LONG Height() const { return m_nHeight; }
	LONG DataWndWidth() const { return m_rectList.Width() - HeaderWidth(); }

	static COLORREF GetNextColor_ThreadPresentationObject();
	static COLORREF GetNextColor_ThreadHeader();

	bool GetTooltipText(CPoint point, CString& strTooltipText) const;
	ThreadsPresentationMode PresentationMode() const { return m_presentationMode; }

	void SetRectangle(CRect rectNew);
	void ZoomInfoChanged();

	void   SetVPosition(int position) { m_vScrollPos = position; }
	void   SetHPosition(double dTimePosition) { m_hScrollPos = dTimePosition; }
	void   GetVScrollInfo(SCROLLINFO& si) const;
	double GetHScrollInfo(SCROLLINFO& si) const;

	double ScrollInfo_Pos() const { return m_hScrollPos; }
	double ScrollInfo_Page() const;
	double ScrollInfo_Max() const;
	double ScrollInfo_PixelRatio() const { return (double)m_rectList.Width() / ScrollInfo_Max(); }
	int    ScrollInfo_TimeToPixel(double dTime) const { return (int)(dTime * ScrollInfo_PixelRatio() + 0.5); }
	double ScrollInfo_PixelToTime(int nPos) const { return (double)nPos / ScrollInfo_PixelRatio(); }

	double GetRelativeTimeFromPixel(CPoint point);
	int GetPixelPositionFromRelativeTime(double time); // returns -1 if time is outside of visible range


	CRect GetDataRect()const;
	CRect GetDrawingRectangle() const;

private:
	void InitColorIndexes();

// Members
protected:
	CRect  m_rectList;
	int    m_vScrollPos;
	double m_hScrollPos;
	std::vector<CThreadPresentationObject> m_threads;
	ThreadsPresentationMode m_presentationMode;

	double m_dStartTime;
	LONG   m_nHeight;
private:
	static size_t m_iColorTPO;
	static size_t m_iColorTH;
};


class CThreadsListFiller
{
private:
	class FillersItem {
		protected:
			std::vector<const DataModel::CThread*> m_threads;

		public:
			FillersItem(const DataModel::CThread* pThread);
			bool InsertItem(const DataModel::CThread* pThread);
			const std::vector<const DataModel::CThread*>& Threads() const { return m_threads; }
	};

protected:
	std::vector<FillersItem> m_items;

public:
	void InsertThread(const DataModel::CThread* pThread);
	void CreateThreadPresentationObjects(LONG& yPosition, LONG nHeaderWidth, CThreadsList* pThreadsList, std::map<CString, COLORREF> & clrFunctions);
};