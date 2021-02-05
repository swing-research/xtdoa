M = 12;
K = 12;
d = 3;
suff = 'reverb_aux';

outDir = strcat('data/',suff);
mkdir(outDir)

N = M + K;  % number of points
J_N = eye(N) - 1/N*ones(N);
solver = "sedumi"; % solver fo cvx

known_delay = [0,0]; % sigma, tau
sdr_init = 1;
rand_restarts = 1; % how many random restarts
c = 343;
rescale =1;

rt60s = [0,0.1,0.2,0.3,0.4,0.5];
R = length(rt60s);
n_runs = 20;

err_lm = zeros(n_runs,length(rt60s)); % final error after refinement
err_tdoa = zeros(n_runs,length(rt60s)); % tdoa estimation error
est_reverb = zeros(n_runs,length(rt60s));

filename = strcat('data/reverb_aux/reverb_M', num2str(M),'_',num2str(d), 'D.mat');
data_all = load(filename);
data_all = data_all.data;

for rti = 1:R
    rt60 = rt60s(rti);
    disp(['rt60 ',num2str(rt60)])  
for ri = 1:n_runs
            disp(['run ',num2str(ri), ' reverb ',num2str(rt60),' (', num2str(M),',' num2str(K),')'])      

           %% Generate points (Forward)
           data_curr = data_all{ri,rti};
           X = data_curr.X; % true positions
           X_c = X * J_N; % centered matrices
           X_r_c = X_c(:, 1:M);
           X_s_c = X_c(:, M+1:end);

           T_true = data_curr.T; % actual TDoA
           
           T = data_curr.T_est; % extracted TDoA
           sqrt(sum((T_true - T).^2, 'all')) / (M*K)

           %% Reconstruct points (Inverse)
           [X_bef, X_aft] = reconstruct_points_rand(T*c, d, known_delay, sdr_init, rand_restarts, solver); % before and after l
         
           %% Compute errors

           % Error after refinement
           X_aft = align_points(X_c, X_aft, rescale); 
           R_aft = X_aft(:, 1:M);
           S_aft = X_aft(:, M+1:end);

           err_lm(ri,rti) = (sum(sqrt(sum((R_aft - X_r_c).^2, 1))) + sum(sqrt(sum((S_aft - X_s_c).^2, 1)))) / N;
        
           err_tdoa(ri,rti) = sum(abs(T_true - T), 'all') / (M*K);
            
           est_reverb(ri,rti) = data_curr.rt60_est;

end
    
end
filename = strcat(outDir,'/matlab_', suff, '_M', num2str(M),'.mat');
save(filename,'err_lm','err_tdoa','est_reverb')