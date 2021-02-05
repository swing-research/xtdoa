Localizing Unsynchronized Sensors with Unknown Sources
================================================================

This repository contains the code and data to reproduce the results of the paper "Localizing Unsynchronized Sensors with Unknown Sources". 
  

Abstract
--------

We propose a method for sensor array self-localization using a set of sources at unknown locations. The sources produce signals whose times of arrival are registered at the sensors. We look at the general case where neither the emission times of the sources nor the reference time frames of the receivers are known. 
Unlike previous work, our method directly recovers the array geometry instead of first estimating the timing information. The key component is a new loss function which is insensitive to the unknown timings. We cast the problem as a minimization of a non-convex functional of the Euclidean distance matrix of microphones and sources subject to certain non-convex constraints. After convexification, we obtain a semidefinite relaxation which gives an approximate solution; subsequent refinement on the proposed loss via the Levenberg-Marquardt scheme gives the final locations. 
Our method achieves state-of-the-art performance in terms of reconstruction accuracy, speed, and ability to work with a small number of sources and receivers. It is also straightforward to add missing measurements, and add arbitrary prior knowledge, for example if either the receiver offsets or the emission times are known, or if the array contains compact subarrays with known geometry.

Authors
-------

Dalia El Badawy, Viktor Larsson, Marc Pollefeys, and Ivan Dokmanić. 


Description
-----------

The `matlab` folder contains the MATLAB scripts and code to reproduce the numerical results. The `python` folder contains the Python scripts to reproduce the plots and figures as well as run the reverberation simulations using Pyroomacoustics.


Dependencies
-------------
### MATLAB

* [MATLAB](https://mathworks.com/)
* [CVX](http://cvxr.com/cvx/download/)


### Python

* [Python 3](https://www.python.org/downloads/)
* [Numpy](http://www.numpy.org/) and [Scipy](https://www.scipy.org/install.html)
* [Pandas](https://pandas.pydata.org)
* [Seaborn](http://seaborn.pydata.org/) and [Matplotlib](http://matplotlib.org/)
* [Pyroomacoustics](https://pypi.org/project/pyroomacoustics/0.1.0/)







