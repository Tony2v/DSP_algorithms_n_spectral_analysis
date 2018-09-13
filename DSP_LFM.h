#pragma once

#include "Parameters.h"
#include "Except.h"

#ifndef DSP_LFM_H
#define DSP_LFM_H
#include <arrayfire.h>
#include <iostream>
#include <string>

using af::Pi;
using std::pow;

class het_exception
{
private:
	SigParams Ref_Sig_Params;
	long long Nin;
	long long Nref;
public:
	het_exception(long long ni, long long nr) : Nin(ni), Nref(nr) {};
	~het_exception() {};

	const std::string what()
	{
		if (Nin != Nref)
			return std::string("Reference signal and input signal dimension mismatch: Nref != Ninput ");
		else return std::string("Smth else hpnd");
	};
};

class Heterodyne
{
private: 
	SigParams Ref_Sig_Params;				// структура с параметрами опорного сигнала
	af::array S_ref;						// generated signal samples

	void SignalGen(const SigParams & sp);	// Throw exceptions on signal type mismatch

public:
	Heterodyne();
	Heterodyne(const SigParams & sp);		// Throw exceptions on signal type mismatch
	~Heterodyne() {};

	void reset(const SigParams & sp);		// Throw exceptions on signal type mismatch
	const SigParams & Heterodyne_State() const { return Ref_Sig_Params; }; 
	
	af::array process_input(const af::array& S_in) const;
};

void Heterodyne::SignalGen(const SigParams & sp)	
{
	// Throw exceptions on signal type mismatch
	// Generates signal if called

	if (sp.Sig == SigType::LFM)
	{
		double N = sp.Trec*sp.Fs;			  // Number of samples

		// Создание массива отсчетов ЛЧМ сигнала без учета сжатия/растяжения из-за Допплера
		double mu = sp.Fdev / sp.Ti;	      // modulation rate parameter for the LFM

		af::array t(af::seq(N) / sp.Fs);	  // creates an array of t [0:N-1]/fdis
		
		af::array Arg_phi = af::array(2.0*Pi*(sp.f0*t + mu*pow(t, 2.0) / 2.0)).as(f64);

		af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
		af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

		S_ref = af::complex(Samp_Re, Samp_Im);  // creates complex LFM signal and copies it to S_ref obj
												// mag of sole real and imag parts become 1/sqrt(2) lesser
												// magnitude of abs(Sref) = 1;
	}
	else
		throw sig_type_mismatch (sp.Sig);
	return;
}


Heterodyne::Heterodyne() 
{
	Ref_Sig_Params.Fdev = 0;
	Ref_Sig_Params.Fs = 0;
	Ref_Sig_Params.Sig = SigType::LFM;
	Ref_Sig_Params.Ti = 0;
	Ref_Sig_Params.Trec = 0;
	S_ref = af::array(0, c64);
}

Heterodyne::Heterodyne(const SigParams & sp)
{
	// Throw exceptions on signal type mismatch
	S_ref = af::array(0, c64);
	if (sp.Sig == SigType::LFM)
	{
		Ref_Sig_Params = sp;
		SignalGen(Ref_Sig_Params);
	}
	else
		throw sig_type_mismatch(sp.Sig);
}

void Heterodyne::reset(const SigParams & sp)
{
	// Throw exceptions on signal type mismatch
	if (sp.Sig == SigType::LFM)
	{
		Ref_Sig_Params = sp;		// reseting signal parameters
		SignalGen(Ref_Sig_Params);
	}
	else
		throw sig_type_mismatch(sp.Sig);
	return;
}

af::array Heterodyne::process_input(const af::array& S_in) const
{		
	if (S_in.elements() != S_ref.elements())
		throw het_exception(S_in.elements(), S_ref.elements());
	
	af::array Sconj = af::conjg(S_ref);  // complex conjugate for the reference signal
	return (S_in * Sconj); // output signal
}

#endif
