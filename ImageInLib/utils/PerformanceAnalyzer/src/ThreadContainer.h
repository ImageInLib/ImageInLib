#pragma once

#include "Thread.h"

namespace DataModel {
	class CThreadContainer
	{
	protected:
		typedef struct Log
		{
			double  time;
			UINT    thread_id;
			UINT    function_id;
			CString log_action;
			CString description;
		};

	public:
		CThreadContainer();
		~CThreadContainer();

		void ProcessLogFile(CString strFileName);
		void Clear(void);

		inline double MinTime() const { return m_dMinTime; }
		inline double MaxTime() const { return m_dMaxTime; }
		inline double MinSize() const { return m_dMinSize; }
		void SetMinTime(double dMinTime) { if (dMinTime < m_dMinTime) m_dMinTime = dMinTime; }
		void SetMaxTime(double dMaxTime) { if (dMaxTime > m_dMaxTime) m_dMaxTime = dMaxTime; }
		void SetMinSize(double dMinRange) { if (dMinRange < m_dMinSize) m_dMinSize = max(5e-8, dMinRange); }

		const std::vector<const CThread*>& Threads() const { return m_threads; }

	protected:
		static bool ParseLog(CString strToParse, Log& logToStore);
		bool AddLog(Log log);

	// members
	protected:
		std::vector<const CThread*> m_threads;
		std::map<UINT, CThread> m_threadsMap;

	protected:
		double m_dMinTime;
		double m_dMaxTime;
		double m_dMinSize;
	};

}