// Security.h: name space Security
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SECURITY_H__924E8C2E_05C4_48D6_8AB5_AECD71BDC145__INCLUDED_)
#define AFX_SECURITY_H__924E8C2E_05C4_48D6_8AB5_AECD71BDC145__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define SECURITY_NAMESPACE_BEGIN  namespace Security {
#define SECURITY_NAMESPACE_END  };

SECURITY_NAMESPACE_BEGIN

class CCertificate;
class CSignature;
class CKey;
class CPublicKey;
class CPrivateKey;
class CSymmetricKey;
class CExchangeKey;
class CKeyRing;
class CBinary;


class CBinary {
public:
	CBinary() : _data(NULL), _size(0) {}
	CBinary(DWORD cbSize) { 
		_data = new BYTE[cbSize]; 
		_size = cbSize; 
	}
	CBinary(BYTE* pData, DWORD cbSize) : _data(NULL), _size(0) { 
		assign(pData, cbSize); 
	}
	CBinary(const CBinary &d) : _data(NULL), _size(0) { 
		assign(d._data, d._size); 
	}
	~CBinary() { 
		if (_data) 
			delete _data; 
	}
	void assign(BYTE* pData, DWORD cbSize) { 
		if (_data) 
			delete _data; 
		_data = new BYTE[cbSize]; 
		memmove(_data, pData, cbSize); 
		_size = cbSize; 
	}
	DWORD size() const { 
		return _size; 
	}
	void resize(DWORD cbSize) { 
		if (cbSize < _size) { 
			_size = cbSize; 
			return; 
		}
		BYTE *p = new BYTE[cbSize];
		memmove(p, _data, _size); 
		_size = cbSize; 
		if (_data) 
			delete _data; 
		_data = p;
	}
	operator const PBYTE() const { 
		return _data; 
	}
	CBinary &operator=(const CBinary &d) { 
		assign(d._data, d._size); 
		return (*this); 
	}
protected:
	BYTE *_data;
	DWORD _size;
};

bool operator<(const Security::CBinary &l, const Security::CBinary &r);
bool operator==(const Security::CBinary &l, const Security::CBinary &r);
bool operator!=(const Security::CBinary &l, const Security::CBinary &r);

SECURITY_NAMESPACE_END



#endif // !defined(AFX_SECURITY_H__924E8C2E_05C4_48D6_8AB5_AECD71BDC145__INCLUDED_)
