#pragma once

#include <map>
#include <vector>
#include "Function.h"

namespace DataModel {
	class CThread
	{
	public:
		CThread(unsigned int nThreadID)
		{
			ThreadID = nThreadID;
			m_dMinTime = FLT_MAX;
			m_dMaxTime = FLT_MIN;
		}

	public:
		std::map<UINT, CFunction> Functions;
		std::vector<UINT> OpenFunctions;
		UINT ThreadID;

		inline double MinTime() const { return m_dMinTime; }
		inline double MaxTime() const { return m_dMaxTime; }
		void UpdateMinTime(double dMinTime) { if (dMinTime < m_dMinTime) m_dMinTime = dMinTime; }
		void UpdateMaxTime(double dMaxTime) { if (dMaxTime > m_dMaxTime) m_dMaxTime = dMaxTime; }

	protected:
		double m_dMinTime;
		double m_dMaxTime;
	};
}
