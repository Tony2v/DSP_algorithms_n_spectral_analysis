// Signal parameters description header file

#pragma once

#include <vector>
#include <valarray>

#ifndef PAR_H
#define PAR_H

#define c_light 299792458
enum class SigType { LFM, PSK };
typedef std::valarray<double> vad;

// структура с параметрами сигнала
struct SigParams
{
	SigType Sig;    // type of the signal (LFM/PSK)
	double Fs;		// sampling rate
	double Ti;		// reflected pulse duration
	double Fdev;	// frequency band
	double f0;		// carrier frequency
	double Trec;	// recieving time
} ;		

enum Mode { BZO, Tracking };	// Processing mode

/*Processing parameters for spectral Analysis Toolset
distances vector - is a 2хN matrix - vector of valarray 2x1 vector object
distance vector must be sorted beforehand in all dimensions from low to higher*/
struct Proc_params
{
	// To be extended ...

	Mode Filt_mode; // Mode: Barrier or Tracking
	std::vector <vad> distances; // low [1] and max [2] bound distances of interval 1-N; 
	double Ti;		// transmitted pulse width
	double Blank;
	double PNI;		// latency where the transmition starts
	double Q;		// duty cycle
	double Fs;		// sampling frequency for processing
	double Norm_factor; 	// normalisation factor for FFT
	double Fband;	// signal's bandwidth

};


#endif