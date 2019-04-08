#ifndef __mutex_h_
#define __mutex_h_

class mutex {

#ifdef _WIN32

	CRITICAL_SECTION c;

public:
    mutex(void) { InitializeCriticalSection(&c); }
    ~mutex(void) { DeleteCriticalSection(&c); }

    void lock(void) const { EnterCriticalSection((LPCRITICAL_SECTION)&c); }
    void unlock(void) const { LeaveCriticalSection((LPCRITICAL_SECTION)&c); }

#else

	pthread_mutex_t     m;
	pthread_mutexattr_t a;

public:

    mutex(void) {
	    pthread_mutexattr_init(&a);
	    pthread_mutexattr_settype(&a, PTHREAD_MUTEX_RECURSIVE_NP);
	    pthread_mutex_init(&m, &a);
    }

    ~mutex(void) {
	    pthread_mutex_destroy(&m);
	    pthread_mutexattr_destroy(&a);
    }

    void lock(void) const { pthread_mutex_lock((pthread_mutex_t*)&m); }
    void unlock(void) const { pthread_mutex_unlock((pthread_mutex_t*)&m); }

#endif

private:
    // dummy copy constructor and operator= to prevent copying
    mutex(const mutex&);
    mutex& operator=(const mutex&);
};


class mutex_lock {
    const mutex& m;
public:
    mutex_lock(mutex& M) : m(M) { m.lock(); }
    mutex_lock(const mutex& M) : m(M) { m.lock(); }
    ~mutex_lock(void) { m.unlock(); }
private:
    // dummy copy constructor and operator= to prevent copying
    mutex_lock(const mutex_lock&);
    mutex_lock& operator=(const mutex_lock&);
};


#endif
