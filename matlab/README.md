Localizing Unsynchronized Sensors with Unknown Sources
================================================================

This folder contains the MATLAB code to reproduce the results of the paper [*Localizing Unsynchronized Sensors with Unknown Sources*](https://arxiv.org/abs/2102.03565).


Description
-----------
### Data

The `data` folder contains the randomly created point configurations for each size stored as `.mat` objects. For example, `data_M7_3D.mat` contains the 200 point configurations for M=K=7 in 3D space. It contains the ToA matrix, the timing offsets, and the true positions.

The real data is stored in `real_data_M12_K65.mat`.

The numerical results are also stored here in subfolders.

### Scripts

One script for each setup at different noise levels and number of microphones. Computes the localization errors and saves the results to a subfolder in the `data` folder.

| Filename | Description|
|-----------|------------|
|`run_exp_lm_rand.m`| our algorithm with or without additional random restarts|
|`run_exp_missing.m`| our algorithm for missing data|
|`run_exp_lm_dist.m`| our algorithms when some distances are known|
|`run_real_exp_K.m`| our algorithm on real data for varying K|
|`run_exp_wang.m`| Wang et al algorithm|
|`run_real_exp_K_wang.m`| Wang et al algorithm on real data for varying K|
|`runtime_comparison.m`| computes the runtime.|
|`run_exp_reverb.m`| our algorithm using TDOA in the presence of reverberation.|
|`run_exp_reverb_wang.m`| Wang et al algorithm using TDOA in the presence of reverberation.|




### Utils

The `utils` folder contains some helper functions.
| Filename | Description|
|-----------|------------|
|`edm.m`| Euclidean Distance Matrix.|
|`align_points.m`| Procrustes. |
|`combine_data.m`| compute input ToA matrix with unknown offsets.|
|`generate_points.m`| creates random points and timing offsets.|
|`generate_dataset.m`| creates and saves 200 configurations for different M=K.|

### Methods

The `methods` folder contains the functions for point reconstruction for different algorithms and variations.

| Filename | Description|
|-----------|------------|
|`relaxcalib_2x2.m`| SDR for solving timing-invariant localization. Can also specify if partially synchronized.|
|`refine_positions_lm.m`| Levenberg-Marquardt (LM) for refinement|
|`relaxcalib_2x2_dist.m`| SDR for solving timing-invariant with some known distances|
|`refine_positions_lm_auglagrange.m`| LM for refinement with some known distances|
|`relaxcalib_2x2_missing.m`| SDR with missing data.|
|`refine_positions_lm_missing.m`| LM with missing data.|
|`estimate_timing_gn.m`| Gauss-Newton algorithm for timing-estimation by Wang et al.|
|`joint_localization_gn.m`| Gauss-Newton algorithm for localization by Wang et al.|
|`reconstruct_points_rand.m`| Calls SDR followed by LM refinement and optionally random restarts|
|`reconstruct_points_dist_augl.m`| Calls SDR followed by LM refinement for known distances|
|`reconstruct_points_missing.m`| Calls SDR followed by LM refinement for missing data|
|`reconstruct_points_wang.m`| Calls Wang's timing estimation then localization.|

