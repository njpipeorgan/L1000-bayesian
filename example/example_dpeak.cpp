#include <random>
#include <iostream>
#include "dpeak.h"

// Inverse Error Function
float myErfInv2(float x) {
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0f : 1.0f;

	x = (1 - x)*(1 + x);        // x = 1 - x*x;
	lnx = logf(x);

	tt1 = 2 / (3.1415*0.147) + 0.5f * lnx;
	tt2 = 1 / (0.147) * lnx;

	return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

int main() {

	//Generate reads


	std::normal_distribution<float> lo_distribution(7.3, 0.27);
	std::normal_distribution<float> hi_distribution(9.4, 0.19);			//The mock distribution is a mixture(2:1) of two normal distributions, centered at 7.3 and 9.4.
	std::uniform_real_distribution<float> bg_distribution(2.0, 14.0);	//The background is obtained by all reads in the well in real data, but we use a uniform distribution for simplicity here. 
	std::discrete_distribution<std::size_t> d({ 32, 64, 4 });			//reads are sampled from a mixture distribution of the peaks and the background with a misidentification rate alpha_c of 4 % .
	std::default_random_engine generator;

	int n_reads = 40;													//sample 40 reads from the mixture distribution

	std::vector<float> values;
	for (int i = 0; i < n_reads; i++)
	{
		switch (d(generator))
		{
		case 0:
			values.push_back(lo_distribution(generator));
			break;
		case 1:
			values.push_back(hi_distribution(generator));
			break;
		case 2:
			values.push_back(bg_distribution(generator));
		}

	}



	// Calculate the probability of being a misidentified bead at the values of all reads. 
	std::vector<float> probbg;
	float pdf_bg = 1. / 12;				// PDF of uniform distribution(2.0,14.0) is a constant 

	for (int i = 0; i < n_reads; i++)
	{
		probbg.push_back(pdf_bg);		// pdf_bg will be replaced by the value of pdf_bg at the value of each read
	}


	// The grid of all possible peak locations for likelihood calculation.
	float grid_start = 2.0;
	float grid_end = 14.0;
	float grid_step = 0.04;

	std::vector<float> peakgrid;
	for (float i = grid_start; i < grid_end; i = i + grid_step) { peakgrid.push_back(i); }
	int gridsize = peakgrid.size();



	float dp52_ratio = 2. / 3;						// ratio of dp52 beads
	float bg = 0.04;								// alpha_c in the paper

	// Calculate the log-likelihood function. 
	std::vector<float> likelihood = dpeak_single(dp52_ratio, bg, values, probbg, peakgrid);
	float max_likelihood = likelihood[0];
	for (auto l : likelihood)
	{
		if (l > max_likelihood)
			max_likelihood = l;
	}
	for (auto& l : likelihood)
	{
		l -= max_likelihood;
	}

	// Marginal likelihood distribution
	std::vector<float> marginal_hi(gridsize, 0.0);
	std::vector<float> marginal_lo(gridsize, 0.0);

	for (int i = 0; i < gridsize; i++)
	{
		for (int j = 0; j < gridsize; j++)
		{
			marginal_hi[i] += std::expf(likelihood[i*gridsize + j]);
			marginal_lo[j] += std::expf(likelihood[i*gridsize + j]);
		}
	}


	// Normalization
	float sum_hi = 0;
	for (auto& n : marginal_hi) sum_hi += n;
	for (auto& n : marginal_hi) n = n / sum_hi;


	float sum_lo = 0;
	for (auto& n : marginal_lo) sum_lo += n;
	for (auto& n : marginal_lo) n = n / sum_lo;


	// Z-score inference

	// In real data, the reference is the distribution of reads from the same gene without perturbation 
	// Here we use a normal distribution (9,1)
	std::vector<float> ref(gridsize, 0.0);
	for (float i = 0; i < gridsize; i = i++)
	{
		float x = 2 + 0.04*i - 9;
		ref[i] = 0.5 * erfc(-x * sqrt(0.5)) - 0.5 * erfc((-x + grid_step) * sqrt(0.5));
	}
	float sum_ref = 0;
	for (auto& n : ref) sum_ref += n;
	for (auto& n : ref) n = n / sum_ref;


	// Calculate z-score for two marginal distributions with ref
	float q1 = 0.0;									
	for (float i = 0; i < gridsize; i = i++)
	{
		float sum_temp = 0;
		for (int j = i; j < gridsize; j++) { sum_temp += marginal_hi[j]; }
		q1 += ref[i] * sum_temp;
	}
	float z1 = sqrt(2) * myErfInv2(2 * q1 - 1);


	float q2 = 0.0;
	for (int i = 0; i < gridsize; i++)
	{
		float sum_temp = 0;
		for (int j = i; j < gridsize; j++) { sum_temp += marginal_lo[j]; }
		q2+= ref[i]*sum_temp;
	}
	float z2 = sqrt(2) * myErfInv2(2 * q2 - 1);

	// z-scores 
	std::cout << "z_hi=" << z1 << "  " << "z_lo=" << z2 << std::endl;
	system("pause");
	return 0;
}
