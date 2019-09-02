// L1000 peak deconvolution based on Bayesian analysis
// 
// Copyright 2019 Tianhuan Lu
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <cmath>

#include <vector>
#include <algorithm>

#include "dpeak.h"


std::vector<float> dpeak_single(
    const float dp52_ratio, 
    const float bg,
    const std::vector<float>& values,
    const std::vector<float>& probbg,
    const std::vector<float>& peakgrid)
{
    const size_t n_read = values.size();
    const size_t n_grid = peakgrid.size();

    if (n_read == 0)
    {
        return std::vector<float>(n_grid * n_grid, 0.0f);
    }

    std::vector<float> probpeak(n_grid * n_read);
    for (size_t g = 0; g < n_grid; ++g)
    {
        float mean = peakgrid[g];
        float dev = PEAK_WIDTH_MODEL(mean);
        for (size_t r = 0; r < n_read; ++r)
        {
            probpeak[g * n_read + r] = 3.30797f * std::pow(3.0f + std::pow((values[r] - mean) / dev, 2), -2) / dev;
            probpeak[g * n_read + r] *= (1.0f - bg);
        }
    }

    std::vector<float> scaled_probbg(probbg);
    for (size_t r = 0; r < n_read; ++r)
        scaled_probbg[r] *= bg;

    std::vector<float> base(n_read);
    std::vector<float> diff(n_read);
    const float inv_num_values = 1.0f / float(n_read);
    const float binomial_var = float(n_read) * dp52_ratio * (1.0f - dp52_ratio);

    std::vector<float> likelihood(n_grid * n_grid);

    for (size_t e52 = 0; e52 < n_grid; ++e52)
    {
        for (size_t e53 = 0; e53 < n_grid; ++e53)
        {
            const float* dp52prob = &probpeak[e52 * n_read];
            const float* dp53prob = &probpeak[e53 * n_read];
            for (size_t r = 0; r < n_read; ++r)
            {
                base[r] = dp52_ratio * dp52prob[r] + (1.0f - dp52_ratio) * dp53prob[r] + scaled_probbg[r];
                diff[r] = inv_num_values * (dp52prob[r] - dp53prob[r]) / base[r];
            }

            float a = 0.0f, b = 0.0f, c = 0.0f;
            for (size_t r = 0; r < n_read; ++r)
            {
                a += std::log(base[r]);
                b += diff[r];
                c += diff[r] * diff[r];
            }
            likelihood[e52 * n_grid + e53] = a + (b * b * binomial_var) / (2.0f + 2.0f * c * binomial_var) 
                                             - 0.5f * std::log(1.0f + c * binomial_var);
        }
    }

    return likelihood;
}

std::vector<float> dpeak_batch(
    const std::vector<float>& dp52_ratio,
    const std::vector<float>& bg,
    const std::vector<std::vector<float>>& values_batch,
    const std::vector<std::vector<float>>& probbg_batch,
    const std::vector<float>& peakgrid)
{
    const size_t grid_size = peakgrid.size();
    const size_t batch_size = values_batch.size();

    std::vector<float> likelihood_batch(batch_size * grid_size * grid_size);

    for (size_t b = 0; b < batch_size; ++b)
    {
        std::vector<float> likelihood = 
            dpeak_single(dp52_ratio[b], bg[b], values_batch[b], probbg_batch[b], peakgrid);
        std::copy_n(likelihood.cbegin(), grid_size * grid_size, 
                    likelihood_batch.begin() + b * grid_size * grid_size);
    }
    return likelihood_batch;
}
