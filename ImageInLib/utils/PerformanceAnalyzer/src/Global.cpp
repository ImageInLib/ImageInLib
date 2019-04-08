#include "stdafx.h"
#include "Global.h"


Global_C *Singleton_C<Global_C>::_pInstance = NULL;
SingletonDestroyer_C<Global_C> Singleton_C<Global_C>::_destroyer;

Global_C::Global_C()
{
}


Global_C::~Global_C()
{
}
