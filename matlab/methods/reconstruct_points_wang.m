function [X_rec] = reconstruct_points_wang(T, dmax, c)
% Reconstruct geometry using Wang et al Algorithm

[delt,eta,X] = estimate_timing_gn(T, dmax);
[R,S,beta0] = joint_localization_gn(T,delt,eta);

X_rec = [R S];
        
