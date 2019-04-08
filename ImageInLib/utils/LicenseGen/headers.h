#ifndef _WIN32
#include <unistd.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#include <pwd.h>
#include <grp.h>
#else
#include <windows.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <signal.h>
#include <string.h>

#ifndef _WIN32
extern "C" {
#include <pthread.h>
extern ssize_t pread __P ((int __fd, __ptr_t __buf, size_t __nbytes,
			   __off_t __offset));
extern ssize_t pwrite __P ((int __fd, __const __ptr_t __buf, size_t __n,
			    __off_t __offset));
extern int pthread_mutexattr_settype __P ((pthread_mutexattr_t *__attr, int __kind));
}
#endif

#ifdef _WIN32
#include <io.h>
#include <process.h>
#include <direct.h>
#endif

# ifdef _MSC_VER
#  pragma warning(disable:4786)
#  pragma warning(disable:4800)
# endif

//STL includes
#include <string>
#include <map>
#include <set>
#include <list>
#include <vector>
#include <algorithm>
#include <stack>
using ::std::for_each;
using ::std::vector;
using ::std::string;
using ::std::list;
using ::std::map;
using ::std::pair;
using ::std::set;
using ::std::stack;


#ifndef _WIN32
typedef CORBA::ULong   DWORD;
typedef CORBA::UShort  WORD;
typedef CORBA::Octet   BYTE;
typedef CORBA::ULong*  PDWORD;
typedef CORBA::UShort* PWORD;
typedef CORBA::Octet*  PBYTE;
typedef void*          HANDLE;
typedef DWORD		   BOOL;
typedef int			   SOCKET;
typedef long long	   LONGLONG;
typedef void		   VOID;
#define closesocket(s) close(s)
#define INVALID_SOCKET -1
#define SOCKET_ERROR   -1
#define TRUE 1
#define FALSE 0
#define MAX_PATH PATH_MAX
#define stricmp strcasecmp
#define strnicmp strncasecmp
#define O_BINARY 0
#else
#define vsnprintf _vsnprintf
#define snprintf _snprintf
#define chdir _chdir
#endif

extern int __argc;
extern char** __argv;

#include "../crypto/crypto.h"
#include "../crypto/auth.h"


