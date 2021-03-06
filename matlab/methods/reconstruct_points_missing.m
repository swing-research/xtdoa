function [X_bef, X_aft] = reconstruct_points_missing(T, d, W, known_delay, sdr_init, rand_restarts, solver)
% Reconstruct geometry when some entries are missing
% Input:
%      T: MxK ToA or TDoA matrix
%      d: dimension of ambient space (2D or 3D)
%      W: mask indicating positions of missing entries
%     known_delay: list indicating if sensors or sources are synchronized
%     [0,0] if fully unsynchronized
%     sdr_init: 1 to use SDR initialization, 0 for random initialization
%     rand_restarts: number of random restarts (1 to only do the LM
%     refinement, 100 to do more)
%     solver: which cvx solver to use ("sedumi")
% Output:
%      X_bef: dx(M+K) point reconstruction after SDR or random init
%      X_aft: dx(M+K) point reconstruction after refinement and optionally
%      restarts
[M,K] = size(T);
if known_delay(1)
    J_M = eye(M);
else  
    J_M = eye(M) - 1/M*ones(M);
end

if known_delay(2)
    J_K = eye(K);
else  
    J_K = eye(K) - 1/K*ones(K);
end

if sdr_init
    % initialize using SDR program
    [X_bef, A_est] = relaxcalib_2x2_missing(T, d, W, known_delay);
else 
    % random initialization

    if d == 2 
        roomdim = [12, 6];
    else
        roomdim = [10, 12, 2.5];
    end

    X_r = diag(roomdim) * rand(d, M);
    X_s = diag(roomdim) * rand(d, K);
    X_bef = [X_r, X_s];
end

% first estimate
R_est = X_bef(:, 1:M);
S_est = X_bef(:, M+1:end);

% refinement
% LM refinement
R_best = R_est; % keep track of points corresponding to lowest residual
S_best = S_est;
err = inf;
dmax = 20^2; % maximum inter-point squared distance in meters sq
for i_rand = 1:rand_restarts
    [R_est, S_est, ~] = refine_positions_lm_missing(R_est, S_est, A_est, T, W, known_delay);

    res = norm(J_M*(sqrt(edm(R_est, S_est)) - T)*J_K, 'fro'); % residual
    X_est = [R_est S_est];
    dmax_est = max(edm(X_est, X_est));
    if res < err & dmax_est<dmax

        err = res;
        R_best = R_est;
        S_best = S_est;
    end


    if err < 1e-6 % change this according to noise level?
        break;
    else
       
        R_est = R_best + 0.3 * randn(size(R_est));
        S_est = S_best + 0.3 * randn(size(S_est));

    end

end

X_aft = [R_best S_best];
            