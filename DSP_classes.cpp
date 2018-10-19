#include "DSP_classes.h"
#include <iostream>

using std::cout;
using std::endl;
using af::Pi;
using std::pow;

Heterodyne::Heterodyne()
{
	Sref_params.Fdev = 0;
	Sref_params.Fs = 0;
	Sref_params.Sig = SigType::LFM;
	Sref_params.Ti = 0;
	Sref_params.Trec = 0;
	Vel = 0;
	cout << "Don't forget to generate the signal by calling SignalGen()!" << endl;

}

Heterodyne::Heterodyne(const SigParams & sp, double Vr)
{
	Vel = Vr;
	Sref_params = sp;
	cout << "Don't forget to generate the signal by calling SignalGen()!" << endl;

}

void Heterodyne::reset(const SigParams & sp, double Vr)
{
	Sref_params = sp;			// resetting signal parameters
	Vel = Vr;					// resetting velocity value
	SignalGen();	
}

af::array Heterodyne::process_input(const af::array& S_in) const
{
	if (S_in.elements() != S_ref.elements())
		throw het_exception(S_in.elements(), S_ref.elements());

	af::array Sconj = af::conjg(S_ref);  // complex conjugate for the reference signal
	return (S_in * Sconj);				 // output signal
}

void Heterodyne::SignalGen()
{
	// Generates signal if called

	double N = Sref_params.Trec*Sref_params.Fs;			  // Number of samples

	af::array t(af::seq(N) / Sref_params.Fs);	  // creates an array of t [0:N-1]/fdis

	double k = (1 - Vel / c_light) / (1 + Vel / c_light); // Doppler's transforming coeffecient

	af::array Arg_phi = af::array(2.0*Pi*Sref_params.f0*k*t).as(f64);

	af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
	af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

	S_ref = af::complex(Samp_Re, Samp_Im);  // creates complex signal and copies it to S_ref obj
											// mag of sole real and imag parts become 1/sqrt(2) lesser
											// magnitude of abs(Sref) = 1;
}

void LFM_Heterodyne::SignalGen()
{
	double N = Sref_params.Trec*Sref_params.Fs;			  // Number of samples

	double mu = Sref_params.Fdev / Sref_params.Ti;	      // modulation rate parameter for the LFM

	af::array t(af::seq(N) / Sref_params.Fs);	  // creates an array of t [0:N-1]/fdis

	double k = (1 - Vel / c_light) / (1 + Vel / c_light); // Doppler's transforming coefficient

	af::array Arg_phi = af::array(2.0*Pi*(Sref_params.f0*k*t + mu*pow(k*t, 2.0) / 2.0)).as(f64);

	af::array Samp_Re = af::cos(Arg_phi);   // generation of the real part of the signal with 1 magnitude!
	af::array Samp_Im = af::sin(Arg_phi);	// generation of the image part of the signal with 1 magnitude!

	S_ref = af::complex(Samp_Re, Samp_Im);  // creates complex LFM signal and copies it to S_ref obj
											// mag of sole real and imag parts become 1/sqrt(2) lesser
											// magnitude of abs(Sref) = 1;
}

void LFM_Heterodyne::reset(const SigParams & sp, double Vr)
{
	Sref_params = sp;		// resetting signal parameters
	Vel = Vr;		// resetting velocity value
	SignalGen();
}
