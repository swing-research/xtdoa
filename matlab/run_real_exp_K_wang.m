clear all;
clc

rng(8);

M = 12;
Ks = 6:11;

n_runs = 200;
d = 3; % ambient dimension
c = 343; % speed of sound

known_delay = [0,0]; % sigma, tau

% suff = 'real_data_sdr_lm';
missing = 0;
minimal = 0;
suff = 'wang_real_data_clean';

rescale = 1; % Procrustes

outDir = strcat('data/',suff);
mkdir(outDir)

dmax = 10;

% load real data
inDir = 'data';
filename = strcat(inDir,'/real_data_M12_K65.mat');
data = load(filename);
T_in = data.T; % available ToAs
[M_tot,K_tot] = size(T_in); % total available sensors and sources
X_r = data.X_r; % mic positions ground truth
W_true = data.W; % mask for outliers

if ~missing
% clean columns only (no outliers)
    s=sum(W_true,1);
    inds = find(s==12);
    T_in = T_in(:,inds);
    [M_tot,K_tot] = size(T_in);
    W_true = W_true(:,inds);
end

for mi=1:length(Ks)
    
    K = Ks(mi);

    N = M + K;  % number of points
    J_N = eye(N) - 1/N*ones(N);
    J_M = eye(M) - 1/M*ones(M);
    
    err_lm = zeros(n_runs, 1); % final error after refinement
    
    for ri = 1:n_runs
        disp(['run ',num2str(ri),' (', num2str(M),',' num2str(K),')'])      

       %% Generate points (Forward)
       % pick random set
       indM = sort(randperm(M_tot,M));
       indK = sort(randperm(K_tot,K));
       
       T = T_in(indM,indK);
       Wm = W_true(indM,indK);
       X_r_curr = X_r(:,indM);
       
       Tdoa = T - T(1,:); % Time-difference-of-arrival
       
     
       %% Reconstruct points (Inverse)
       if ~missing
        X_aft = reconstruct_points_wang(Tdoa, dmax, c);
       end
       %% Compute errors

       
       % Error after refinement
       R_aft = X_aft(:, 1:M);
       S_aft = X_aft(:, M+1:end);
       R_aft = align_points(X_r_curr, R_aft, rescale); % only align microphones    
       

       err_lm(ri) = sum(sqrt(sum((R_aft - X_r_curr).^2, 1))) / M;


    end
    % save results
    filename = strcat(outDir,'/matlab_', suff, '_M', num2str(M), '_K', num2str(K), '.mat');
    save(filename,'err_lm')
end
  