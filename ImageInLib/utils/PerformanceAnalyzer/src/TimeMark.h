#pragma once

#include "stdafx.h"
#include "zoominfo.h"

#define TIME_NOT_SET -1.0

class TimeMark
{
public:
	TimeMark(DataModel::CZoomInfo& zoomInfo, COLORREF color) 
		: m_zoomInfo(zoomInfo)
		, m_color(color)
		, m_dTime(TIME_NOT_SET)
	{ }

	TimeMark(const TimeMark & t)
		: m_zoomInfo(t.m_zoomInfo)
		, m_color(t.m_color)
		, m_dTime(t.m_dTime)
	{ }

	~TimeMark() { }

	TimeMark & operator=(const TimeMark & t)
	{
		m_zoomInfo = t.m_zoomInfo;
		m_color = t.m_color;
		m_dTime = t.m_dTime;
	}

	enum Style {
		Normal,
		LeftArrow,
		RightArrow
	};

	void Draw(CDC* pDC, CRect rect, Style style = Normal);

	bool SetTime(double dTime);
	double GetTime() const { return m_dTime; }
	bool IsValid() const { return m_dTime != TIME_NOT_SET; }
	void Clear() { m_dTime = TIME_NOT_SET; }
	
protected:
	double m_dTime;
	DataModel::CZoomInfo& m_zoomInfo;
	COLORREF m_color;
};


