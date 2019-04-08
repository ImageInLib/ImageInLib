#if !defined(__HEADERS_H__)
#define __HEADERS_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#pragma warning(disable:4786)
#pragma warning(disable:4800)

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#include <windows.h>
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <stdarg.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <limits.h>
#include <direct.h>

//STL includes
#pragma warning(push, 3)
#include <string>
using ::std::string;
using ::std::wstring;
#include <exception>
using ::std::exception;
#include <map>
using ::std::map;
using ::std::pair;
using ::std::multimap;
#include <set>
using ::std::set;
#include <list>
using ::std::list;
#include <vector>
using ::std::vector;
#include <stack>
using ::std::stack;
#include <queue>
using ::std::queue;
#include <algorithm>
#include <limits>
#pragma warning(pop)

#define __APPNAME__ "rcdiff_Hotkey"
#define __APP_VERSION__ "0.002"

#endif // !defined(__HEADERS_H__)
