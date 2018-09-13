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
private:
	SigGenPar prm;

public:
	Sig_Generator(const SigGenPar& sp) : prm(sp) {};
	~Sig_Generator() {};

	const af::array Get(const double& Tr) const;

};

const af::array Sig_Generator::Get(const double& Tr) const
{
	if (prm.Sig == SigType::LFM)
	{
		double mu = prm.Fdev / prm.Ti;	// modulation rate parameter for the LFM

		if (Tr < prm.Ti + prm.Blank)
		{
			double T_input = Tr - prm.Blank;			// rec signal duration
			double N = floor(T_input*prm.Fs);			// number of signal samples
			af::array t(af::seq(N) / prm.Fs);			// creates an array of t [0:N-1]/fdis
			double T_shift = prm.Ti + prm.Blank - Tr;	// time shift for the rec signal phase

			// Phase of the rec signal:
			af::array Arg_phi = af::array(2.0*Pi*(prm.f0*(t + T_shift) + mu*pow((t + T_shift), 2.0) / 2.0)).as(f64);

			af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
			af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

			af::array S = af::complex(Samp_Re, Samp_Im);	//  creates complex LFM signal

			// Padding the input sequence with End_pad zeros to get the real input:
			double End_pad = ceil(((prm.Q - 1)*prm.Ti - T_input)*prm.Fs);
			
			af::array ZZend = af::constant(0., End_pad, c64);

			S = af::join(0, S, ZZend);

			return S;
		}
		else if (Tr < (prm.Q - 1)*prm.Ti + prm.Blank)
		{
			double N = floor(prm.Ti*prm.Fs);				// number of signal samples

			af::array t(af::seq(N) / prm.Fs);	// creates an array of t [0:N-1]/fdis

												// Phase of the input signal:
			af::array Arg_phi = af::array(2.0*Pi*(prm.f0*t + mu*pow(t, 2.0) / 2.0)).as(f64);

			af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
			af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

			af::array S = af::complex(Samp_Re, Samp_Im);	//  creates complex LFM signal

			// Padding the signal with zeros at the beginning and the end:
			double Beg_pad = floor((Tr - prm.Ti - prm.Blank)*prm.Fs);
			double End_pad = floor(((prm.Q - 1)*prm.Ti + prm.Blank - Tr)*prm.Fs);

			af::array ZZbeg = af::constant(0., Beg_pad, c64);
			af::array ZZend = af::constant(0., End_pad, c64);

			S = af::join(0, ZZbeg, S);
			S = af::join(0, S, ZZend);

			return S;
		}
		else
		{
			double T_input = prm.Q*prm.Ti + prm.Blank - Tr;	// rec signal duration
			double N = floor(T_input*prm.Fs);				// number of signal samples
			af::array t(af::seq(N) / prm.Fs);	// creates an array of t [0:N-1]/fdis

												// Phase of the input signal:
			af::array Arg_phi = af::array(2.0*Pi*(prm.f0*t + mu*pow(t, 2.0) / 2.0)).as(f64);

			af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
			af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

			af::array S = af::complex(Samp_Re, Samp_Im);	//  creates complex LFM signal

			double Beg_pad = floor((Tr - prm.Ti - prm.Blank)*prm.Fs);
			af::array ZZbeg = af::constant(0., Beg_pad, c64);
			S = af::join(0, ZZbeg, S);

			return S;
		}
	}
	else if (prm.Sig == SigType::PSK)
	{
		// Returning void af::array object - has not been developed yet
		return af::array();
	}
	else std::cerr << "Invalid type of the signal!\n";
}

#endif