#pragma once

#include "arrayfire.h"
#include "Parameters.h"
#include <iostream>

#ifndef SIG_GEN_H
#define SIG_GEN_H

struct SigGenPar
{
	SigType Sig;	// type of the signal (LFM/PSK)
	double Ti;		// transmitted pulse width
	double Blank;	
	double Q;		// duty cycle
	double Fs;		// sampling frequency for processing
	double f0;		// carrier frequency
	double Fdev;	// bandwidth
};

using af::Pi;

class Sig_Generator
{
public:
	Sig_Generator(const SigGenPar& sp) : prm(sp) {};
	virtual ~Sig_Generator() {};

	virtual const af::array Get(double Tr, double Vr);

protected:
	void get_the_strobe(double Tr, double k);

	SigGenPar prm;	
	double Beg_pad;		// number of zeros to pad the beginning
	double End_pad;		// number of zeros to pad the end
	double T_shift;		// time shifting of the phase for the signal to be generated
	double T_input;		// rec signal duration
	
	af::array Arg_phi;	// phase of the signal to be generated with paramters prm and Tr and k = (1-Vr/c)/(1+Vr/c)
};

void Sig_Generator::get_the_strobe(double Tr, double k)
{
	if (Tr < prm.Ti + prm.Blank)
	{
		T_input = k*(Tr - prm.Blank);		// rec signal duration

		T_shift = prm.Ti + prm.Blank - Tr;	// time shift for the rec signal phase		

		// Padding the input sequence with End_pad zeros to get the real input:
		End_pad = ceil(((prm.Q - 1)*prm.Ti - T_input)*prm.Fs);
	}
	else if (Tr < (prm.Q - k)*prm.Ti + prm.Blank)
	{
		T_input = k*prm.Ti;			// time duration of the received input signal

		T_shift = 0;

		// Padding the signal with zeros at the beginning and the end:
		Beg_pad = floor((Tr - prm.Ti - prm.Blank)*prm.Fs);
		End_pad = floor((prm.Q*prm.Ti + prm.Blank - Tr - T_input)*prm.Fs);
	}
	else
	{
		T_input = prm.Q*prm.Ti + prm.Blank - Tr;	// rec signal duration
		T_shift = 0;
		Beg_pad = floor((Tr - prm.Ti - prm.Blank)*prm.Fs);
	}
}


const af::array Sig_Generator::Get(double Tr, double Vr)
{	
	double k = (1 - Vr / c_light) / (1 + Vr / c_light); // коэффициент трасформации временного масштаба сигнала из-за Допплера
	
	get_the_strobe(Tr, k);
	
	af::array ZZbeg = af::constant(0., Beg_pad, c64);
	af::array ZZend = af::constant(0., End_pad, c64);
	
	double N = floor(T_input*prm.Fs);			// number of signal samples
	af::array t(af::seq(N) / prm.Fs);			// creates an array of t [0:N-1]/fdis
	
	switch (prm.Sig)
	{
	case SigType::LFM:
	{
		double mu = prm.Fdev / prm.Ti;	// modulation rate parameter for the LFM

		// Phase of the rec signal:
		Arg_phi = af::array(2.0*Pi*(prm.f0*(k*t + T_shift) + mu*pow((k*t + T_shift), 2.0) / 2.0)).as(f64);
	}
		break;
	
	case SigType::PSK:
		// Returning void af::array object - has not been developed yet
		return af::array();
		break;

	default:
		std::cerr << "Invalid type of the signal!\n";
		break;
	}
		
	af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
	af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

	af::array S = af::complex(Samp_Re, Samp_Im);	//  creates complex LFM signal

	S = af::join(0, ZZbeg, S);
	S = af::join(0, S, ZZend);

	return S;
}

#endif