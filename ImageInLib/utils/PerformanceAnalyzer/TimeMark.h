#pragma once

#include "stdafx.h"
#include "zoominfo.h"

class TimeMark
{
public:
	TimeMark(DataModel::CZoomInfo& m_zoomInfo);
	~TimeMark(){};

	void Draw(CDC* pDC, CRect rect);
	bool SetTime(double dTime);
	double GetTime()const{ return m_dTime; }
	
protected:
	double m_dTime;
	DataModel::CZoomInfo& m_zoomInfo;
};


