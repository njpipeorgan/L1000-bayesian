
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>

#include "../../src/dpeak.h"


// the inverse error function
float erf_inv(float x);

int main()
{
    // Step 1: Generate the log2-FI values for a bead color.

    // The mock distributions for two peaks centered at 7.3 and 9.4. 
    // They will be mixed at a ratio of 2:1 (dp52:dp53).
    std::normal_distribution<float> hi_distribution(7.3, 0.27);
    std::normal_distribution<float> lo_distribution(9.4, 0.19);

    // The background distribution should be obtained from all reads in the 
    // well, but we use a uniform distribution for simplicity here. 
    std::uniform_real_distribution<float> bg_distribution(2.0, 14.0);

    // The reads are sampled from a mixture distribution of the peaks and the 
    // background with a misidentification rate alpha_c of 4%.
    std::discrete_distribution<> mixture({ 0.64, 0.32, 0.04 });

    // Sample 40 reads from the mixture distribution
    size_t n_reads = 40;
    std::default_random_engine generator;

    std::vector<float> values;
    for (size_t i = 0; i < n_reads; ++i)
    {
        switch (mixture(generator))
        {
        case 0: values.push_back(hi_distribution(generator)); break;
        case 1: values.push_back(lo_distribution(generator)); break;
        case 2: values.push_back(bg_distribution(generator));
        }
    }


    // Step 2: Get all other parameters for likelihood calculation.

    // Calculate the probability of coming from a misidentified bead at the 
    // values of all reads. In this case, the PDF is a constant: 1/12. 
    std::vector<float> prob_bg(n_reads, 1.0f / 12);

    // The grid of all possible peak locations for likelihood calculation: 
    // equally spaced between 2.0 and 14.0 with a step of 0.03
    float grid_start = 2.0;
    float grid_step = 0.03;
    float grid_size = 400;

    std::vector<float> peak_grid;
    for (size_t i = 0; i < grid_size; ++i)
        peak_grid.push_back(grid_start + float(i) * grid_step);

    // the ratio of dp52 beads
    float dp52_ratio = 2. / 3;

    // the ratio of background beads. In the actual calculations, multiple 
    // ratios are used, and the likelihood functions are combined by their 
    // probability, but we only use a single ratio (= alpha_c) for simplicity. 
    float bg_ratio = 0.04;


    // Step 3: Calculate the log-likelihood function.

    // Call the peak deconvolution function.
    std::vector<float> likelihood = dpeak_single(dp52_ratio, bg_ratio, values, prob_bg, peak_grid);

    // Adjust the log-likelihoods to improve the accuracy of 
    // marginal distributions later.
    float max_likelihood = *std::max_element(likelihood.cbegin(), likelihood.cend());
    for (float& l : likelihood)
        l -= max_likelihood;

    
    // Step 4: Calculate the marginal distributions. 

    // Calculate the marginal distributions.
    std::vector<float> marginal_hi(grid_size, 0.0);
    std::vector<float> marginal_lo(grid_size, 0.0);

    for (int i = 0; i < grid_size; ++i)
    {
        for (int j = 0; j < grid_size; ++j)
        {
            float p = std::exp(likelihood[i * grid_size + j]);
            marginal_hi[i] += p;
            marginal_lo[j] += p;
        }
    }

    // Normalize the marginal distributions
    float sum_hi = std::accumulate(marginal_hi.cbegin(), marginal_hi.cend(), 0.f);
    for (auto& m : marginal_hi)
        m /= sum_hi;

    float sum_lo = std::accumulate(marginal_lo.cbegin(), marginal_lo.cend(), 0.f);
    for (auto& m : marginal_lo)
        m /= sum_lo;


    // Step 5: z-score inference. 

    // In real data, the reference is the averaged marginal distribution of 
    // the same gene in all wells from the same plate (i.e. plate control). 
    // Here, we employ an idealized reference distribution N(mu=9.0, var=1.0). 
    auto ref_cdf = [](float x) { return 0.5 * std::erfc((9.0f - x) / std::sqrt(2.0f)); };
    std::vector<float> ref;
    for (size_t i = 0; i < grid_size; ++i)
    {
        float x = (2.0f + grid_step * float(i));
        ref.push_back(ref_cdf(x) - ref_cdf(x - grid_step));
    }

    // Calculate quantile for the high-abundance peak against the reference 
    // distribution above. 
    float q_hi = 0.0;									
    for (size_t i = 0; i < grid_size; ++i)
    {
        float prob_above = std::accumulate(marginal_hi.cbegin() + i, marginal_hi.cend(), 0.0f);
        q_hi += ref[i] * prob_above;
    }

    // Calculate the z-score according to the quantile. 
    float zscore_hi = std::sqrt(2.0f) * erf_inv(2.0f * q_hi - 1.0f);

    std::cout << "zscore_hi = " << zscore_hi << "\n";

    return 0;
}


float erf_inv(float x)
{
    float sgn = (x < 0) ? -1.0f : 1.0f;
    float lnx = std::log(1.0f - x * x);
    float tt1 = 2 / (3.1415 * 0.147) + 0.5f * lnx;
    float tt2 = 1 / (0.147) * lnx;
    return sgn * std::sqrt(-tt1 + std::sqrt(tt1 * tt1 - tt2));
}
