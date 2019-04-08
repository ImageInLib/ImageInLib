#include "stdafx.h"
#include "ThreadContainer.h"
#include "TextFileReader.h"

using namespace DataModel;

CThreadContainer::CThreadContainer()
{
	Clear();
}


CThreadContainer::~CThreadContainer()
{
}

//reads log file (with format defined in CPerformanceLog) and fill container Threads
void CThreadContainer::ProcessLogFile(CString strFileName)
{
	CTextFileReader fileReader(strFileName);
	if (!fileReader.OpenFile()) {
		AfxMessageBox(_T("Error occured while opening log file."));
		return;
	}

	CString str;
	for (int nResult = 0; nResult == 0;)
	{
		// read line from log file
		nResult = fileReader.ReadLine(str);
		// if nResult == 0, line was read; if nResult == -1, end of line was reached; else an error occured
		if (nResult != 0) {
			if (nResult > 0) {
				AfxMessageBox(_T("Error occured while reading log file."));
				m_threads.clear();
			}
			break;
		}
		// parse it
		Log log;
		if (!ParseLog(str, log))
		{
			AfxMessageBox(_T("Incorrect file format."));
			m_threads.clear();
			break;
		}
		// add data from from log
		if (!AddLog(log))
		{
			AfxMessageBox(_T("Incorrect file format."));
			m_threads.clear();
			break;
		}
	};

	fileReader.CloseFile();
}

void CThreadContainer::Clear(void)
{
	m_threads.clear();
	m_threadsMap.clear();
	m_dMinTime = FLT_MAX;
	m_dMaxTime = FLT_MIN;
	m_dMinSize = FLT_MAX;
}

bool CThreadContainer::ParseLog(CString strToParse, Log& logToStore)
{
	if (strToParse.IsEmpty())
		return false;

	Log log;
	CString str;
	int iStart = 0;
	int iEnd = 0;
	
	// find time stamp
	iEnd = strToParse.Find(TCHAR('\t'), iStart);
	if (iEnd == -1)
		return false;
	errno = 0;
	str = strToParse.Mid(iStart, iEnd);
	log.time = _ttof(str);
	if ((log.time == 0 && errno != 0) || log.time < 0)
		return false;
	// find thread ID
	iStart = iEnd + 1;
	iEnd = strToParse.Find(TCHAR('\t'), iStart);
	if (iEnd == -1)
		return false;
	errno = 0;
	str = strToParse.Mid(iStart, iEnd - iStart);
	log.thread_id = (UINT)_ttol(str);
	if (log.thread_id == 0 && errno != 0)
		return false;
	// find function ID
	iStart = iEnd + 1;
	iEnd = strToParse.Find(TCHAR('\t'), iStart);
	if (iEnd == -1)
		return false;
	errno = 0;
	str = strToParse.Mid(iStart, iEnd - iStart);
	log.function_id = (UINT)_ttol(str);
	if (log.function_id == 0 && errno != 0)
		return false;
	// find action string
	iStart = iEnd + 1;
	iEnd = strToParse.Find(TCHAR('\t'), iStart);
	if (iEnd == -1)
		return false;
	log.log_action = strToParse.Mid(iStart, iEnd - iStart);
	if (log.log_action.IsEmpty())
		return false;
	// the rest is description
	iStart = iEnd + 1;
	iEnd = strToParse.Find(TCHAR('\t'), iStart);
	if (iEnd == -1) {
		iEnd = strToParse.GetLength();
		if (iStart > iEnd)
			return false;
	}
	log.description = strToParse.Mid(iStart, iEnd - iStart);

	logToStore = log;
	return true;
}

bool CThreadContainer::AddLog(Log log)
{
	bool bAddNewFunction = true;
	const UINT& nThreadID = log.thread_id;
	const UINT& nFunctionID = log.function_id;
	std::map<UINT, CFunction>::iterator itFunction;
	std::map<UINT, CThread>::iterator itThread;

	// if thread specified by nThreadID doesn't exist, create it
	itThread = m_threadsMap.find(nThreadID);
	if (itThread == m_threadsMap.end())
	{
		itThread = m_threadsMap.insert(std::pair<UINT, CThread>(nThreadID, CThread(nThreadID))).first;
		m_threads.push_back(&itThread->second);
	}
	else
	{
		itFunction = itThread->second.Functions.find(nFunctionID);
		if (itFunction != itThread->second.Functions.end())
			bAddNewFunction = false;
	}
	// if function doesn't exist, create it
	if (bAddNewFunction)
		itFunction = itThread->second.Functions.insert(std::pair<UINT, CFunction>(nFunctionID, CFunction(nThreadID))).first;
	
	// set data to function object
	if (log.log_action == "start")
	{
		itFunction->second.Start = log.time;
		itFunction->second.Level = itThread->second.OpenFunctions.size();
		itFunction->second.Description = log.description;
		itThread->second.UpdateMinTime(log.time);
		SetMinTime(log.time);
		itThread->second.OpenFunctions.push_back(nFunctionID);
	}
	else if (log.log_action == "end")
	{
		auto & openFunctions = itThread->second.OpenFunctions;
		if (openFunctions.empty() || openFunctions.back() != nFunctionID)
		{
			// function end without function start or incorrect nesting of functions
			ASSERT(0);
			return false;
		}
		openFunctions.pop_back();
		itFunction->second.Stop = log.time;
		itThread->second.UpdateMaxTime(log.time);
		SetMaxTime(log.time);
		SetMinSize(itFunction->second.Stop - itFunction->second.Start);
	}
	else
	{
		ASSERT(0); //NotImplemented;
		return false;
	}
	return true;
}
