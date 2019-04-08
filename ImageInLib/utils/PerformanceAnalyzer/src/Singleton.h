#pragma once

/**********************************************************
 * Singleton classes declaration & definition             *
 **********************************************************/

template <class SingletonType>
class SingletonDestroyer_C
{
private:
	SingletonType *_pDoomed;
public:
	SingletonDestroyer_C(SingletonType *pDoomed = NULL) {
		_pDoomed = pDoomed;
	};
	virtual ~SingletonDestroyer_C(void) {
		if (_pDoomed)
			delete _pDoomed;
	};
	void SetDoomed(SingletonType *pDoomed) {
		_pDoomed = pDoomed;
	};
};


template <class Type>
class Singleton_C
{
	friend SingletonDestroyer_C < Type > ;
private:
	static Type *_pInstance;
	static SingletonDestroyer_C<Type> _destroyer;
protected:
	Singleton_C(void) {}
	virtual ~Singleton_C(void) {}
public:
	static Type *Instance(void) {
		if (!_pInstance) {
			_pInstance = new Type;
			_destroyer.SetDoomed(_pInstance);
		}
		return _pInstance;
	};
	static Type *GetCurrentInstance() { return _pInstance; }
	static void DestroyInstance(void) {
		if (_pInstance) {
			delete _pInstance;
			_pInstance = NULL;
			_destroyer.SetDoomed(NULL);
		}
	};
};