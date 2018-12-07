// Defines Spectral Analysis toolset for the signal processing after 
// heterodyning the input
// !!! check the Parametres structure before issuing !!!
// exceptions has not yet been established

#include "Parameters.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <arrayfire.h>
#include <iostream>
#include <valarray>
#include <tuple>
#include <utility>

#pragma once
#ifndef SA_ALGO_H
#define SA_ALGO_H

using std::vector;
using af::Pi;

// Functor to derive the size of FFT for input data based on input time value
class nfft_pow2
{
private:
	double Fs;
public:
	nfft_pow2(const double& f) : Fs(f) {};
	~nfft_pow2() {};
	double operator()(const std::pair<double,double>& t) const { return pow(2, ceil(log2(2*(t.second-t.first)*Fs))); };
};

// Function to obtain latency values of input signal from distance value
template <typename T>
T dist_2_latency(const T& d) { return d*2.0 / c_light; };
template <> vad dist_2_latency(const vad& d) { return d*2.0 * (1.0 / c_light); };

double input_t_out_of_latency(const vad& lat) { return (lat[1] - lat[0]); };

// Main class for Spectral Analysis toolset
// To be developed as well..
// exception-free for a while
class SA
{

private:
	Proc_params SA_settings;									// Structure with a parameters' set for SA
	vector<af::array> Freq_grid;								// âåêòîð, ñîäåðæàùèé ñåòêó ÷àñòîò â êàæäîé èç ïîëîñ ñïåêòðîàíàëèçàòîðà
	vector<std::tuple<double, double, double>> F_min_max_step;	// ìèí/ìàêñ ÷àñòîòà â êàæäîé èç ïîëîñ è øàã ïî ÷àñòîòå

public:
	vector<af::array> output_SA; // âûõîä ÑÀ

	const vector<af::array>& get_freq_grid() const { return Freq_grid; };
	const vector<std::tuple<double, double, double>>& get_F_min_max_step() const { return F_min_max_step; };
	
	SA(const Proc_params& pr) : SA_settings(pr) {};  // constructor of SA
	~SA() {};

	const vector <af::array>& process(const af::array& in);
	
};

const vector <af::array>&  SA::process(const af::array& in)
{
	unsigned long dist_size = SA_settings.distances.size();	// number of the distance slots
		
	// Function to obtain pair<double, double> input start/stop time based on latency values derived from serving distances vector
	auto input_from_latency = [&](const vad& lat)
	{
		double t_begin;
		double t_end;
		if (lat[0] <= SA_settings.Ti + SA_settings.Blank)
		{
			t_begin = 0.;
			t_end = (lat[1] < (SA_settings.Q - 1)*SA_settings.Ti + SA_settings.Blank) ? (lat[1] - SA_settings.Blank) : (SA_settings.Q - 1)*SA_settings.Ti;
			return std::make_pair(t_begin, t_end);
		}
		else
		{
			t_begin = (lat[0] - SA_settings.Ti - SA_settings.Blank);
			t_end = (lat[1] < (SA_settings.Q - 1)*SA_settings.Ti + SA_settings.Blank) ? (lat[1] - SA_settings.Blank) : (SA_settings.Q - 1)*SA_settings.Ti;
		}
		return std::make_pair(t_begin, t_end);
	};

	if (SA_settings.Filt_mode == Mode::BZO)
	{
		double mu = SA_settings.Fband / SA_settings.Ti; // ñêîðîñòü ìîäóëÿöèè äëÿ Ë×Ì ñèãíàëà

		vector<vad>::iterator it_b = SA_settings.distances.begin();
		vector<vad>::iterator it_e = SA_settings.distances.end();

		vector<double> nfft_size(dist_size);	// fft size to filter that data input accordingly with corresponding distance slot
		vector<vad> latency(dist_size);			// latency of received signal corresponding to distance values
		vector<std::pair<double,double>> t_input(dist_size);		// input time slots (just time lengths) accordingly with corresponding distance values

		std::transform(it_b, it_e, latency.begin(), dist_2_latency<vad>);  //obtaining latency values of input signal from distance values
		std::transform(latency.begin(), latency.end(), t_input.begin(), input_from_latency); // obtaining time slots (just time duration period) from latency max and min values
		std::transform(t_input.begin(), t_input.end(), nfft_size.begin(), nfft_pow2(SA_settings.Fs)); // nfft size vector corresponding to time slots to filter received data

		double N_strobe;	// size of input strobe;
					
		for (int i = 0; i < dist_size; i++)
		{
			double begin = ceil(t_input[i].first*SA_settings.Fs);
			double end = floor(t_input[i].second*SA_settings.Fs) - 1;

			// Writing input strobe's sample sequence:
			af::seq input_samples(begin, end);			// sequence's numbers of input samples
			af::array input_strobe = in(input_samples);		// acquiring the samples of the strobe from input signal

			// Weighting with Hamming window:
			N_strobe = end - begin + 1;

			af::array n = af::array(af::seq(N_strobe)).as(f64);
			af::array Hamming_win = (0.54 - 0.46*af::cos(2.0 * Pi * n / (N_strobe - 1))).as(f64);

			//dim_t dH = Hamming_win.elements();
			//dim_t dI = input_strobe.elements();

			input_strobe = input_strobe * Hamming_win;		// Now weighted

			af::array input_spectre = af::fftNorm(input_strobe, SA_settings.Norm_factor, nfft_size[i]);	// returns input_strobe's spectre values in number nfft_size[i] after FFT with normalisation factor;
			
			double Fmin_high = (SA_settings.Ti - (latency[i][0] - SA_settings.Blank))*mu; // ÷àñòîòà ìèíèìëüíîé îáðàáàòûâàåìîé äàëüíîñòè äëÿ äàííîãî ñòðîáà ïðèåìà (âûøå ÷åì ÷àñòîòà ìàêñ äàëüíîñòè)
			double Fmax_low = (SA_settings.Ti - (latency[i][1] - SA_settings.Blank))*mu;  // ÷àñòîòà ìàêñ îáðàáàòûâàåìîé äàëüíîñòè äëÿ äàííîãî ñòðîáà ïðèåìà (íèæå ÷åì ÷àñòîòà ìèí äàëüíîñòè)
			
			double k = nfft_size[i] / N_strobe; // êîýôôèöèåíò, ó÷èòûâàþùèé âî ñêîëüêî ðàç äîïîëíåííàÿ íóëÿìè ïîñë-òü îòñ÷åòîâ áîëüøå ÷åì êîë-âî ïîñòóïàþùèõ îò÷åòîâ â ñòðîáå ïðèåìà

			double Fstep = 1 / (t_input[i].second - t_input[i].first) / k; // øàã ïî ñåòêå ÷àñòîò ÑÀ

			F_min_max_step.push_back(std::make_tuple(Fmin_high, Fmax_low, Fstep));

			double min_R_sample = Fmin_high * (t_input[i].second - t_input[i].first) * k;  // íîìåð îòñ÷åòà ÷àñòè ñïåêòðà Ñ êîòîðîãî òðåáóåòñÿ ïðîèçâîäèòü çàïèñü, äâèãàÿñü â ñòîðîíó íóëÿ
			double max_R_sample = Fmax_low * (t_input[i].second - t_input[i].first) * k;   // íîìåð ïîñëåäíåãî îòñ÷åòà ÷àñòè ñïåêòðà ÄÎ êîòîðîãî òðåáóåòñÿ ïðîèçâîäèòü çàïèñü, äâèãàÿñü â ñòîðîíó íóëÿ
			
			af::seq ind_to_write = af::seq(min_R_sample, max_R_sample, -1);

			Freq_grid.push_back(af::array(ind_to_write).as(f64)*Fstep);

			output_SA.push_back(input_spectre(ind_to_write)); // reading only those values which needed based on high and low strobe' frequency values 
						
			//double mm = af::max<double>(af::abs(input_spectre));

		}

		return output_SA;
	}
	else
	{
		std::cout << "Tracking has not yet been developed!" << std::endl;
		return output_SA;
	}
		
};

#endif

/*// Functor to obtain input time values from latency of the signal based on distance vector, Ti, Blank and duty cycle Q
class input_t_out_of_latency
{
private:
double Ti;
double Blank;
double Q;

public:
input_t_out_of_latency(const double t, const double b, const double q) : Ti(t), Blank(b), Q(q) {};
~input_t_out_of_latency() {};

double operator()(const double lat)
{
std::initializer_list<double> ilist = { lat - Blank , Q*Ti + Blank - lat, Ti };
return std::min<double>(ilist);
};
};*/
