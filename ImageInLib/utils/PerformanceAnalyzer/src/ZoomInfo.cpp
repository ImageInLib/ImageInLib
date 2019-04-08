#include "stdafx.h"
#include "ZoomInfo.h"

using namespace DataModel;

CZoomInfo::CZoomInfo()
{
	Reset();
}

void CZoomInfo::Reset()
{
	m_dOriginalZoomCoef = 1.;
	m_dActualZoomCoef = 1.;
	m_dZoomRangeStart = 0.;
	m_dZoomRangeWidth = _100ns;
	m_dZoomRangeMinWidth = _100ns;
	m_dDataWidth = 0.;
	m_dDataShift = 0.;
}

CZoomInfo::~CZoomInfo()
{
}

bool CZoomInfo::SetDataShift(double dDataShift)
{
	if (dDataShift < 0)
		dDataShift = 0;
	if (dDataShift == m_dDataShift)
		return false;
	m_dDataShift = dDataShift;
	return true;
}

void CZoomInfo::UpdateZoomCoefs(double dDataWndWidth, double dDataWidth)
{
	if (dDataWndWidth < 0)
		dDataWndWidth = 0;
	m_dDataWidth = dDataWidth;
	m_dOriginalZoomCoef = (dDataWidth != 0) ? dDataWndWidth / dDataWidth : 1.;
	UpdateActualZoom();
}

double CZoomInfo::SetZoomRange(double dRangeStart, double dRangeWidth)
{
	// new range should be within data range
	if (dRangeWidth > m_dDataWidth)
		dRangeWidth = m_dDataWidth;
	if (dRangeWidth < m_dZoomRangeMinWidth)
		dRangeWidth = m_dZoomRangeMinWidth;
	if ((dRangeStart + dRangeWidth) > m_dDataWidth)
		dRangeStart = m_dDataWidth - dRangeWidth;
	if (dRangeStart < 0)
		dRangeStart = 0;
	// is there something to change?
	if (m_dZoomRangeStart == dRangeStart && m_dZoomRangeWidth == dRangeWidth)
		return m_dZoomRangeStart;
	// set zoom ramge
	m_dZoomRangeStart = dRangeStart;
	m_dZoomRangeWidth = dRangeWidth;
	// now calculate actual zoom coef.
	UpdateActualZoom();
	return m_dZoomRangeStart;
}

double CZoomInfo::SetMinZoomRange(double dMinRange)
{
	if (dMinRange < _100ns)
		dMinRange = _100ns;
	if (m_dZoomRangeMinWidth != dMinRange) {
		m_dZoomRangeMinWidth = dMinRange;
		if (m_dZoomRangeWidth > m_dZoomRangeMinWidth)
			SetZoomRange(m_dZoomRangeStart, m_dZoomRangeMinWidth);
	}
	return m_dZoomRangeMinWidth;
}

bool CZoomInfo::UpdateActualZoom()
{
	double dZoom = (m_dZoomRangeWidth > 0) ? m_dDataWidth * m_dOriginalZoomCoef / m_dZoomRangeWidth : m_dOriginalZoomCoef;
	if (dZoom == m_dActualZoomCoef)
		return false;
	m_dActualZoomCoef = dZoom;
	return true;
}
