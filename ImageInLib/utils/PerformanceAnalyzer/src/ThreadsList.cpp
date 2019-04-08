#include "stdafx.h"
#include "ThreadsList.h"
#include "Global.h"
#include <algorithm>


size_t CThreadsList::m_iColorTPO = 0;
size_t CThreadsList::m_iColorTH = 0;

CThreadsList::CThreadsList() :
	m_dStartTime(0.),
	m_presentationMode(TPM_Separately),
	m_vScrollPos(0),
	m_hScrollPos(0.),
	m_nHeight(0)
{
	InitColorIndexes();
}


CThreadsList::~CThreadsList()
{
}

bool CThreadsList::Initialize(ThreadsPresentationMode presentationMode)
{
	InitColorIndexes();
	m_presentationMode = presentationMode;
	m_nHeight = 0;
	m_vScrollPos = 0;
	m_hScrollPos = 0.;

	// at first, obtain data range (min time, max time) and update zoom info
	m_dStartTime = FLT_MAX;
	double dMaxTime = 0;
	const std::vector<const DataModel::CThread*>& threads = Global_C::Instance()->ThreadContainer().Threads();
	for (const DataModel::CThread* thread : threads)
	{
		if (thread->MinTime() < m_dStartTime)
			m_dStartTime = thread->MinTime();
		if (thread->MaxTime() > dMaxTime)
			dMaxTime = thread->MaxTime();
	}
	if (m_dStartTime == FLT_MAX)
		m_dStartTime = 0;
	DataModel::CZoomInfo& zoomInfo = Global_C::Instance()->ZoomInfo();
	double dDataRange = dMaxTime - m_dStartTime;
	zoomInfo.SetMinZoomRange(Global_C::Instance()->ThreadContainer().MinSize() * 2.);
	zoomInfo.SetDataShift(m_dStartTime);
	zoomInfo.UpdateZoomCoefs(DataWndWidth(), dDataRange);
	zoomInfo.SetZoomRange(0, dDataRange);
	zoomInfo.SetDataShift(m_dStartTime);

	// display each thread as separate thread presentation object
	std::map<CString, COLORREF> clrFunctions;
	if (m_presentationMode == TPM_Separately)
	{
		for (const auto& thread : threads)
		{
			CThreadsListFiller filler;
			filler.InsertThread(thread);
			filler.CreateThreadPresentationObjects(m_nHeight, HeaderWidth(), this, clrFunctions);
		}
	}
	else // display non-overlapping threads jointly
	{
		CThreadsListFiller filler;
		for (const auto& thread : threads)
		{
			filler.InsertThread(thread);
		}
		filler.CreateThreadPresentationObjects(m_nHeight, HeaderWidth(), this, clrFunctions);
	}
	return true;
}

void CThreadsList::Clear()
{
	m_threads.clear();
	Global_C::Instance()->ZoomInfo().Reset();
}

void CThreadsList::SetRectangle(CRect rectNew)
{
	if (m_rectList == rectNew)
		return;
	m_rectList = rectNew;
	// set zoom
	DataModel::CZoomInfo& zoomInfo = Global_C::Instance()->ZoomInfo();
	zoomInfo.UpdateZoomCoefs(DataWndWidth(), zoomInfo.DataWidth());
	ZoomInfoChanged();
}

CRect CThreadsList::GetDataRect()const
{
	CRect rc = m_rectList;
	rc.left += HeaderWidth();
	rc.bottom = rc.top + m_nHeight;
	return rc;
}

CRect CThreadsList::GetDrawingRectangle() const
{
	return CRect(m_rectList.TopLeft(), CSize(m_rectList.Width(), Height()));
}

void CThreadsList::ZoomInfoChanged()
{
	for (CThreadPresentationObject& thread : m_threads)
	{
		thread.Invalidate();
	}
}

LONG CThreadsList::HeaderWidth() const
{
	return (m_presentationMode == TPM_Separately) ? 70 : 0;
}

double CThreadsList::ScrollInfo_Page() const
{
	return Global_C::Instance()->ZoomInfo().ZoomRange();
}

double CThreadsList::ScrollInfo_Max() const
{
	return Global_C::Instance()->ZoomInfo().DataWidth();
}

void CThreadsList::GetVScrollInfo(SCROLLINFO& si) const
{
	int nVertScrollPos = m_vScrollPos;

	si.cbSize = sizeof(SCROLLINFO);
	si.fMask = SIF_DISABLENOSCROLL | SIF_POS | SIF_PAGE | SIF_RANGE;
	si.nMin = 0;
	si.nMax = Height();
	si.nPage = m_rectList.Height() + 1;
	si.nPos = (nVertScrollPos + (int)si.nPage > si.nMax) ? max(0, si.nMax - (int)si.nPage + 1) : nVertScrollPos;
}

double CThreadsList::GetHScrollInfo(SCROLLINFO &si) const
{
	double dPos = m_hScrollPos;
	double dPage = ScrollInfo_Page();
	double dMax = ScrollInfo_Max();
	if (dPos + dPage > dMax)
		dPos = max(0., dMax - dPage);
	si.cbSize = sizeof(SCROLLINFO);
	si.fMask = SIF_DISABLENOSCROLL | SIF_POS | SIF_PAGE | SIF_RANGE;
	si.nMin = 0;
	si.nMax = m_rectList.Width();
	si.nPage = (UINT)ScrollInfo_TimeToPixel(dPage) + 1;
	si.nPos = ScrollInfo_TimeToPixel(dPos);
	return dPos;
}


void CThreadsList::Draw(CDC &dc)
{
	CRgn rgn;
	rgn.CreateRectRgnIndirect(m_rectList);
	dc.SelectClipRgn(&rgn, RGN_AND);

	CSize size(m_rectList.Width(),m_rectList.Height());
	for (auto& thread : m_threads) {
		if (thread.Draw(dc, size, m_rectList.TopLeft(), m_hScrollPos, m_vScrollPos) > 0)
			break;
		dc.SelectClipRgn(&rgn, RGN_OR);
	}
}

void CThreadsList::InitColorIndexes()
{
	m_iColorTPO = m_iColorTH = 0;
	CFunctionPresentationObject::InitColorIndexes();
}

COLORREF CThreadsList::GetNextColor_ThreadPresentationObject()
{
	const COLORREF colors[] = { RGB(245, 245, 245), RGB(250, 250, 250) };

	m_iColorTPO++;
	if (m_iColorTPO >= sizeof(colors) / sizeof(colors[0]))
		m_iColorTPO = 0;
	return colors[m_iColorTPO];
}

COLORREF CThreadsList::GetNextColor_ThreadHeader()
{
	const COLORREF colors[] = { RGB(215, 215, 215), RGB(230, 230, 230) };

	m_iColorTH++;
	if (m_iColorTH >= sizeof(colors) / sizeof(colors[0]))
		m_iColorTH = 0;
	return colors[m_iColorTH];
}

double CThreadsList::GetRelativeTimeFromPixel(CPoint point)
{
	int xOffset = point.x - m_rectList.left - HeaderWidth();
	return  Global_C::Instance()->ZoomInfo().PixelPositionToTimePosition(xOffset) + m_hScrollPos;
}

int CThreadsList::GetPixelPositionFromRelativeTime(double time)
{
	double dTimeOffset = time - m_hScrollPos;
	double pixelPosition = m_rectList.left + HeaderWidth() + floor(Global_C::Instance()->ZoomInfo().TimePositionToPixelPosition(dTimeOffset) + 0.5);
	if (pixelPosition < HeaderWidth() || pixelPosition > m_rectList.right)
		return -1;
	return (int)pixelPosition;
}

bool CThreadsList::GetTooltipText(CPoint point, CString& strTooltipText) const
{
	if (!m_rectList.PtInRect(point))
		return false;

	point -= m_rectList.TopLeft();

	for (const CThreadPresentationObject& thread : m_threads)
	{
		if (thread.GetTooltipText(point, strTooltipText, m_vScrollPos))
			break;
	}

	return true;
}

CThreadsListFiller::FillersItem::FillersItem(const DataModel::CThread* pThread)
{
	m_threads.push_back(pThread);
}

bool CThreadsListFiller::FillersItem::InsertItem(const DataModel::CThread* pThread)
{
	double dMin = pThread->MinTime();
	double dMax = pThread->MaxTime();
	std::vector<const DataModel::CThread*>::iterator itWhereToInsert;
	for (itWhereToInsert = m_threads.begin(); itWhereToInsert != m_threads.end(); itWhereToInsert++)
	{
		if (dMax < (*itWhereToInsert)->MinTime())
			break;
		else if (dMin > (*itWhereToInsert)->MaxTime())
			continue;
		else
			return false;
	}
	m_threads.insert(itWhereToInsert, pThread);
	return true;
}

void CThreadsListFiller::InsertThread(const DataModel::CThread* pThread)
{
	std::vector<FillersItem>::iterator it = m_items.begin();
	for (; it != m_items.end(); it++)
	{
		if (it->InsertItem(pThread))
			return;
	}
	m_items.push_back(FillersItem(pThread));
}

void CThreadsListFiller::CreateThreadPresentationObjects(LONG& yPosition, LONG nHeaderWidth, CThreadsList* pThreadsList, std::map<CString, COLORREF> & clrFunctions)
{
	for (auto& item : m_items)
	{
		const std::vector<const DataModel::CThread*>& threads = item.Threads();
		// format thread header string if there is only one thread in the line
		CString strHeader;
		if (threads.size() == 1)
			strHeader.Format(_T("%d"), threads[0]->ThreadID);
		// create thread presentation object
		CThreadPresentationObject threadPresentationObject(
			yPosition,
			nHeaderWidth,
			pThreadsList->CThreadsList::GetNextColor_ThreadHeader(),
			strHeader,
			pThreadsList->CThreadsList::GetNextColor_ThreadPresentationObject());
		// add functions
		for (auto pThread : threads)
		{
			// collect functions
			std::vector<const DataModel::CFunction*> functions;
			functions.reserve(pThread->Functions.size());
			for (const auto& function : pThread->Functions)
			{
				if (function.second.Start <= function.second.Stop)
					functions.push_back(&function.second);
			}
			// add them to the thread presentation object
			threadPresentationObject.AddFunctions(clrFunctions, functions);
		}
		threadPresentationObject.Invalidate();
		pThreadsList->m_threads.push_back(threadPresentationObject);
		yPosition += threadPresentationObject.Height();
	}
}
