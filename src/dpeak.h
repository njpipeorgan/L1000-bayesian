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

// scale parameter of the peak as a function of expression level x
#define PEAK_WIDTH_MODEL(x) (0.15f + 5.3f * std::exp2f(-0.75f * x))


// Calculate the log-likelihood function of the position of the peaks. 
// The result (gridsize x gridsize) is flattened into one-dimension.

std::vector<float> dpeak_single(
// the ratio of dp52 bead set
	const float dp52_ratio,
// the ratio of reads from color misidentification
	const float bg,
// log2-FI values as the reads
	const std::vector<float>& values,
// the PDF of background (all reads in the same well) at each read
	const std::vector<float>& probbg,
// the locations of the peaks for likelihood calculation
	const std::vector<float>& peakgrid);


// Calculate for multiple beads in batch. 
// The result is of size (N x gridsize x gridsize) and flattened. 

std::vector<float> dpeak_batch(
	const std::vector<float>& dp52_ratio,
	const std::vector<float>& bg,
	const std::vector<std::vector<float>>& values_batch,
	const std::vector<std::vector<float>>& probbg_batch,
	const std::vector<float>& peakgrid);
	