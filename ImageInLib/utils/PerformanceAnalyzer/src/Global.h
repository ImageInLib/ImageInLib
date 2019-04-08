#pragma once

#include "Singleton.h"
#include "ThreadContainer.h"
#include "ZoomInfo.h"

class Global_C : public Singleton_C<Global_C>
{
	friend Singleton_C<Global_C>;
	friend SingletonDestroyer_C<Global_C>;

public:
	Global_C();
	virtual ~Global_C();

protected:
	DataModel::CThreadContainer m_threadContainer;
	DataModel::CZoomInfo        m_zoomInfo;

public:
	DataModel::CThreadContainer& ThreadContainer() { return m_threadContainer; }
	DataModel::CZoomInfo&        ZoomInfo()        { return m_zoomInfo; }
};

