Localizing Unsynchronized Sensors with Unknown Sources
================================================================

This folder contains the Python code to plot the figures in the paper as well as run the reverberation simulations using Pyroomacoustics.

Description
-----------
### Scripts

Each script loads the appropriate results from the `matlab/data` folder and plots a figure stored in `figures`.


| Filename | Description|
|-----------|------------|
|`plot_compare.py`| Fig. 2, 4.|
|`plot_runtime.py`| Fig. 3. Runtime comparison between our approach and Wang et al's.|
|`plot_missing.py`| Fig. 5. Results for our algorithm with missing data.|
|`simulate_room.py`| Pyroomacoustics simulation and TDOA estimation.|
|`plot_reverb.py`| Fig. 6. Results in the presence of reverberation.|
|`plot_real_data_input.py`| Fig. 7. Illustration of real data. |
|`plot_real_data.py`|  Fig. 8. Real data results for varying K.|
|`plot_real_data_missing.py`| Fig. 9. Real data results with missing entries. |
|`aux_tdoa.py`| Code by Robin Scheibler. TDOA estimation. |






