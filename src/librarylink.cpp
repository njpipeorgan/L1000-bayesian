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

#include <algorithm>
#include <vector>

#include "WolframLibrary.h"


std::vector<float> dpeak_batch(
    const std::vector<float>&, const std::vector<float>&, const std::vector<std::vector<float>>& values_batch,
    const std::vector<std::vector<float>>&, const std::vector<float>& peakgrid);

EXTERN_C DLLEXPORT int WolframLibrary_initialize(WolframLibraryData lib_data)
{
    return 0;
}

EXTERN_C DLLEXPORT void WolframLibrary_uninitialize(WolframLibraryData)
{
}

EXTERN_C DLLEXPORT int wll_dpeak(WolframLibraryData libdata, mint argc, MArgument* args, MArgument res)
{
    MTensor mt_values    = MArgument_getMTensor(args[0]);
    MTensor mt_dp52ratio = MArgument_getMTensor(args[1]);
    MTensor mt_bgratio   = MArgument_getMTensor(args[2]);
    MTensor mt_probbg    = MArgument_getMTensor(args[3]);
    MTensor mt_peakgrid  = MArgument_getMTensor(args[4]);

    const mreal* values_ptr    = libdata->MTensor_getRealData(mt_values);
    const mreal* dp52ratio_ptr = libdata->MTensor_getRealData(mt_dp52ratio);
    const mreal* bgratio_ptr   = libdata->MTensor_getRealData(mt_bgratio);
    const mreal* probbg_ptr    = libdata->MTensor_getRealData(mt_probbg);
    const mreal* peakgrid_ptr  = libdata->MTensor_getRealData(mt_peakgrid);

    size_t batch_size = libdata->MTensor_getDimensions(mt_values)[0];
    size_t values_width = libdata->MTensor_getDimensions(mt_values)[1];
    size_t grid_size = libdata->MTensor_getDimensions(mt_peakgrid)[0];

    if (batch_size   != libdata->MTensor_getDimensions(mt_dp52ratio)[0] ||
        batch_size   != libdata->MTensor_getDimensions(mt_bgratio)[0]   ||
        batch_size   != libdata->MTensor_getDimensions(mt_probbg)[0]    ||
        values_width != libdata->MTensor_getDimensions(mt_probbg)[1])
    {
        return LIBRARY_DIMENSION_ERROR;
    }

    std::vector<std::vector<float>> values;
    std::vector<std::vector<float>> probbg;
    std::vector<float> dp52ratio(dp52ratio_ptr, dp52ratio_ptr + batch_size);
    std::vector<float> bgratio(bgratio_ptr, bgratio_ptr + batch_size);
    std::vector<float> peakgrid(peakgrid_ptr, peakgrid_ptr + grid_size);

    for (size_t b = 0; b < batch_size; ++b)
    {
        const mreal* values_begin = values_ptr + b * values_width;
        const mreal* probbg_begin = probbg_ptr + b * values_width;
        size_t num_values = std::find_if(values_begin, values_begin + values_width, [](mreal v) { return v < 0.0; }) - values_begin;

        values.push_back(std::vector<float>(values_begin, values_begin + num_values));
        probbg.push_back(std::vector<float>(probbg_begin, probbg_begin + num_values));
    }

    std::vector<float> likelihood = dpeak_batch(dp52ratio, bgratio, values, probbg, peakgrid);

    mint likelihood_dims[3] ={ mint(batch_size), mint(grid_size), mint(grid_size) };
    MTensor mt_likelihood;
    libdata->MTensor_new(MType_Real, 3, likelihood_dims, &mt_likelihood);
    mreal* likelihood_ptr = libdata->MTensor_getRealData(mt_likelihood);

    std::copy_n(likelihood.cbegin(), likelihood.size(), likelihood_ptr);
    MArgument_setMTensor(res, mt_likelihood);

    return LIBRARY_NO_ERROR;
}

