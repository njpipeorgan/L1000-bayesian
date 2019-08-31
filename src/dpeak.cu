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

#include <vector>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// maximum number of reads per bead color
#define CUDA_PROBPEAK_STRIDE 512

// scale parameter of the peak as a function of expression level x
#define PEAK_WIDTH_MODEL(x) (0.15f + 5.3f * exp2f(-0.75f * x))


__global__ void probpeak_kernel(const float bg1m, const float* values, const float* peakgrid, float* probpeak)
{
    int g = blockIdx.x;
    int i = threadIdx.x;

    float mean = peakgrid[g];
    float dev = PEAK_WIDTH_MODEL(mean);
    float inv_dev = 1.0f / dev;
    float prob = (values[i] - mean) * inv_dev;
    probpeak[g * CUDA_PROBPEAK_STRIDE + i] = bg1m * 3.30797f * powf(3.0f + prob * prob, -2.0f) * inv_dev;
}

__global__ void prob_kernel(const unsigned int num_values, unsigned int n_reduction,
    const float inv_num_values, const unsigned int grid_size,
    const float dp52_ratio, const float dp53_ratio,
    const float* probpeak, const float* probbg, float* intparams)
{
    int e52 = blockIdx.x;
    int e53 = blockIdx.y;
    int r = threadIdx.x;

    const float* dp52prob = probpeak + (e52 * CUDA_PROBPEAK_STRIDE);
    const float* dp53prob = probpeak + (e53 * CUDA_PROBPEAK_STRIDE);
    float* totalparams = intparams + (e52 * grid_size + e53) * 3;

    extern __shared__ float params[];

    float base = dp52_ratio * dp52prob[r] + dp53_ratio * dp53prob[r] + probbg[r];
    float diff = inv_num_values * (dp52prob[r] - dp53prob[r]) / base;

    params[r * 3 + 0] = __logf(base);
    params[r * 3 + 1] = diff;
    params[r * 3 + 2] = diff * diff;
    __syncthreads();

    if (r < int(num_values - n_reduction))
    {
        params[r * 3 + 0] += params[(r + n_reduction) * 3 + 0];
        params[r * 3 + 1] += params[(r + n_reduction) * 3 + 1];
        params[r * 3 + 2] += params[(r + n_reduction) * 3 + 2];
    }
    __syncthreads();
    n_reduction >>= 1;
    for (; n_reduction > 0; n_reduction >>= 1)
    {
        if (r < n_reduction)
        {
            params[r * 3 + 0] += params[(r + n_reduction) * 3 + 0];
            params[r * 3 + 1] += params[(r + n_reduction) * 3 + 1];
            params[r * 3 + 2] += params[(r + n_reduction) * 3 + 2];
        }
        __syncthreads();
    }
    if (r == 0)
    {
        totalparams[0] = params[0];
        totalparams[1] = params[1];
        totalparams[2] = params[2];
    }
}

__global__ void reduce_kernel(const float binomial_var, const float* intparams, float* likelihood)
{
    int e52 = blockIdx.x;
    int e53 = threadIdx.x;
    int stride = blockDim.x;

    const float* params_ptr = intparams + (e52 * stride + e53) * 3;
    float* likelihood_ptr = likelihood + (e52 * stride + e53);

    float a = params_ptr[0];
    float b = params_ptr[1];
    float c = params_ptr[2];
    *likelihood_ptr = a + (b * b * binomial_var) / (2.0f + 2.0f * c * binomial_var) - 0.5f * __logf(1.0f + c * binomial_var);
}

void cudadpeak_single(
    const float dp52_ratio, const float bg,
    const std::vector<float>& values, float* d_values,
    float* h_probpeak, float* d_probpeak,
    std::vector<float> probbg, float* h_probbg, float* d_probbg,
    const std::vector<float>& peakgrid, float* d_peakgrid,
    float* d_intparams, float* h_likelihood, float* d_likelihood)
{
    const size_t num_values = values.size();
    const size_t grid_size = peakgrid.size();

    if (num_values == 0)
    {
        std::fill_n(h_likelihood, grid_size * grid_size, 0.0f);
        return;
    }

    for (size_t i = 0; i < num_values; ++i)
        h_probbg[i] = probbg[i] * bg;
    cudaMemcpy(d_probbg, h_probbg, sizeof(float) * num_values, cudaMemcpyHostToDevice);
    cudaMemcpy(d_values, values.data(), sizeof(float) * num_values, cudaMemcpyHostToDevice);

    probpeak_kernel<<<grid_size, num_values>>>(1.0f - bg, d_values, d_peakgrid, d_probpeak);

    float inv_num_values = 1.0f / float(num_values);
    float binomial_var = float(num_values) * dp52_ratio * (1.0f - dp52_ratio);
    unsigned int n_reduction = 1;
    while (n_reduction * 2 < num_values)
        n_reduction *= 2;

    dim3 dim_grid(grid_size, grid_size);
    dim3 dim_block(num_values);
    size_t shared_memory_size = sizeof(float) * num_values * 3;
    prob_kernel<<<dim_grid, dim_block, shared_memory_size>>>(
        int(num_values), n_reduction, inv_num_values, int(grid_size),
        dp52_ratio, 1.0f - dp52_ratio, d_probpeak, d_probbg, d_intparams);

    std::vector<float> h_intparams(10000);
    cudaMemcpy(h_intparams.data(), d_intparams, h_intparams.size() * 4, cudaMemcpyDeviceToHost);

    reduce_kernel<<<grid_size, grid_size>>>(binomial_var, d_intparams, d_likelihood);
    cudaMemcpy(h_likelihood, d_likelihood, sizeof(float) * grid_size * grid_size, cudaMemcpyDeviceToHost);
}

std::vector<float> cudadpeak_batch(
    const std::vector<float>& dp52_ratio,                // ratio of dp52 beads
    const std::vector<float>& bg,                        // alpha_c in the paper
    const std::vector<std::vector<float>>& values_batch, // log2-FI in batches
    const std::vector<std::vector<float>>& probbg_batch, // PDF of all reads at each log2-FI
    const std::vector<float>& peakgrid)                  // grid points of peak locations
{
    const size_t batch_size = values_batch.size();
    const size_t grid_size = peakgrid.size();

    std::vector<float> probpeak(grid_size * CUDA_PROBPEAK_STRIDE);
    std::vector<float> probbg(CUDA_PROBPEAK_STRIDE);
    std::vector<float> likelihood_batch(batch_size * grid_size * grid_size);

    float* h_probpeak = probpeak.data();
    float* h_probbg = probbg.data();

    float* d_probpeak;
    float* d_values;
    float* d_probbg;
    float* d_peakgrid;
    float* d_intparams;
    float* d_likelihood;

    cudaMalloc(&d_probpeak,   sizeof(float) * grid_size * CUDA_PROBPEAK_STRIDE);
    cudaMalloc(&d_values,     sizeof(float) * CUDA_PROBPEAK_STRIDE);
    cudaMalloc(&d_probbg,     sizeof(float) * CUDA_PROBPEAK_STRIDE);
    cudaMalloc(&d_peakgrid,   sizeof(float) * grid_size);
    cudaMalloc(&d_intparams,  sizeof(float) * grid_size * grid_size * 3);
    cudaMalloc(&d_likelihood, sizeof(float) * grid_size * grid_size);

    cudaMemcpy(d_peakgrid, peakgrid.data(), sizeof(float) * grid_size, cudaMemcpyHostToDevice);

    for (size_t b = 0; b < batch_size; ++b)
    {
        float* h_likelihood = likelihood_batch.data() + b * grid_size * grid_size;
        cudadpeak_single(dp52_ratio[b], bg[b],
            values_batch[b], d_values,
            h_probpeak, d_probpeak,
            probbg_batch[b], h_probbg, d_probbg,
            peakgrid, d_peakgrid,
            d_intparams, h_likelihood, d_likelihood);
    }

    cudaFree(d_likelihood);
    cudaFree(d_intparams);
    cudaFree(d_peakgrid);
    cudaFree(d_probbg);
    cudaFree(d_values);
    cudaFree(d_probpeak);

    return likelihood_batch;
}

