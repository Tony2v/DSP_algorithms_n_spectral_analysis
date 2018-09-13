#pragma once

#ifndef EXCEPT_H
#define EXCEPT_H

#include "Parameters.h"
#include <iostream>

class sig_type_mismatch
{
private:
	char Sig_names[2][4] = { "LFM", "PSK" };
public:
	SigType signame;	// if the exception is thrown signal type Flag is copied here

	sig_type_mismatch(SigType s) : signame(s) {};

	inline void mesg()
	{
		std::cout << "Heterodine is not defined for the "
			<< Sig_names[static_cast <int> (signame)] << std::endl;
	};
};


#endif