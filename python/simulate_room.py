#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pyroomacoustics as pra
from scipy.spatial import distance
from matplotlib import pyplot as plt
from scipy.io import savemat, loadmat
from scipy.io.wavfile import read
from pathlib import Path
import os
from aux_tdoa import aux_tdoa, parabola
import itertools

###### Simulate different reverberation times and estimate TDOAs

np.random.seed(88)


M = 12 # number of microphones
K = 12 # number of sources
N = M + K # total number of points
d = 3 # ambient dimension
fs = 48000 # sampling rate
c = 343 # speed of sound

rt60s = [0, 0.1,0.2, 0.3, 0.4, 0.5] # desired reverberation times in seconds
runs = 20 # per reverb time

snr = 15 # simulation SNR


path = Path(__file__).parent / os.path.join('..','matlab','data') # path to the saved results from matlab
outpath = os.path.join(path,'reverb_aux')
if not os.path.exists(outpath):
    os.makedirs(outpath)


# random receiver and source positions
roomdim = np.array([10, 10, 3])  

# set params for each reverbertation time 3D
rtparam = {}
rtparam[0] = [1.,0]
rtparam[0.1] = [0.7,15]
rtparam[0.2] = [0.5,23]
rtparam[0.3] = [0.37,30]
rtparam[0.4] = [0.3,35]
rtparam[0.5] = [0.25,37]

data = {}
truth = {}
thelist = []
interp = 2
phat = False
for ri in range(runs):
    # generate points
    X_r = np.diag(roomdim).dot(np.random.rand(d, M))
    X_s = np.diag(roomdim).dot(np.random.rand(d, K))
    X = np.hstack((X_r, X_s))

    # generate K non-overlapping source signals
    T = 0.5 # duration of sound within slot in seconds
    T_tot = T + 3 # duration of slot in seconds
    
    N = np.int32(T*fs) # number of samples
    N_tot = np.int32(T_tot*fs)
    offsets = np.zeros((K,)) # ground truth emission times
    signals = np.zeros((K,N))
    for ki in range(K):
        noise = np.random.randn(N,) # signal of duration T
        signals[ki,:] = noise
        offset = np.random.uniform(0,2) # start time within first two seconds
        offsets[ki] = offset

    # True TDoA matrix
    D_true = distance.cdist(X_r.T,X_s.T,'euclidean')
    T = distance.cdist(X_r.T,X_s.T,'euclidean')/c
    T += (np.ones((M, 1)).dot(offsets[:,np.newaxis].T))
    T = T - T[0,:] # TDoA wrt first mic
    
    revlist = []
    for ti,t60 in enumerate(rt60s):
        abso,max_ord = rtparam[t60]
    
        # create room
        room = pra.ShoeBox(roomdim, fs=fs, max_order=max_ord, absorption=abso) # specify reverb here

        mic = pra.MicrophoneArray(X_r, fs)
        room.add_microphone_array(mic)

        # generate K non-overlapping source signals
        for ki in range(K):
            offset = offsets[ki]
            tot_offset = ki*T_tot + offset
            room.add_source(X_s[:,ki], signal=signals[ki,:], delay=tot_offset)

        room.image_source_model()
        room.simulate(snr=snr)

        avg_rt60 = np.mean(room.measure_rt60())
        print(t60,'estimated rt60',avg_rt60)
        
        # Extract TDoA
        T_est = np.zeros((M,K))
        for mi in range(M):
            for ki in range(K):
                # extract kth source segment in mth microphone
                seg_ref = room.mic_array.signals[0,ki*N_tot:(ki+1)*N_tot] # zeroth mic is reference
                seg = room.mic_array.signals[mi,ki*N_tot:(ki+1)*N_tot]
    
                n_iter = 10
                f_init_kwargs = {"phat": phat, "interp": interp, "fs": 1.0}
                tau_init = parabola(seg, seg_ref, **f_init_kwargs)
                tau0 = aux_tdoa(seg, seg_ref, phat=phat, t0=tau_init, n_iter=n_iter,)

                T_est[mi,ki] = tau0 / fs
                
        err_tdoa = np.mean(np.abs(T-T_est))
        print(ri,ti,t60,'tdoa err',err_tdoa)

        curr = {}
        curr['T'] = T
        curr['X'] = X
        curr['offsets'] = offsets
        curr['T_est'] = T_est
        curr['rt60'] = t60
        curr['rt60_est'] = avg_rt60
        
        revlist.append(curr)
        
    thelist.append(revlist)

data_mat = {'data':thelist}
filename = os.path.join(outpath,'reverb_M%s_%sD.mat'%(M,d))
savemat(filename, data_mat)
