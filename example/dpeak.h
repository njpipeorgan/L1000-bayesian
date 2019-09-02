
#include <cmath>

#include <vector>
#include <algorithm>


std::vector<float> dpeak_single(
	const float dp52_ratio, const float bg,
	const std::vector<float>& values,
	const std::vector<float>& probbg,
	const std::vector<float>& peakgrid);

std::vector<float> dpeak_batch(
	const std::vector<float>& dp52_ratio,                // ratio of dp52 beads
	const std::vector<float>& bg,                        // alpha_c in the paper
	const std::vector<std::vector<float>>& values_batch, // log2-FI in batches
	const std::vector<std::vector<float>>& probbg_batch, // PDF of all reads at each log2-FI
	const std::vector<float>& peakgrid);                  // grid points of peak locations