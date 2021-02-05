clear all;
clc

% Compare runtimes

rng(8);

Ms = 7:12;
n_runs = 20;
d = 3; % ambient dimension
c = 343; % speed of sound
dmax = 10;

suff = 'runtime_comp_err';

known_delay = [0,0]; % sigma, tau
stdev = 1e-5;

outDir = strcat('data/',suff);
mkdir(outDir)


rescale = 1;

runtime_ours = zeros(length(Ms),n_runs);
runtime_wang = zeros(length(Ms),n_runs);

err_ours = zeros(length(Ms),n_runs);
err_wang = zeros(length(Ms),n_runs);
for mi=1:length(Ms)
    
    M = Ms(mi); % number of microphones
    K = M; % number of sources
    filename = strcat('data/data_M', num2str(M),'_',num2str(d), 'D.mat');
    data = load(filename);
    data = data.data;
    
    N = M + K;  % number of points
    J_N = eye(N) - 1/N*ones(N);   
    
    for ri = 1:n_runs
        disp(['run ',num2str(ri), ' noise ',num2str(stdev),' (', num2str(M),',' num2str(K),')'])
        
        %% Generate points (Forward)
        X = data{ri,4}; % true positions
        X_c = X * J_N; % centered matrices
        X_r_c = X_c(:, 1:M);
        X_s_c = X_c(:, M+1:end);
        
        T = combine_data(data{ri,1}, data{ri,2}, data{ri,3}, stdev, c, known_delay);
        
        %% Reconstruct points (Inverse)
        Tdoa = T - T(1,:);

        tic
        X_aft = reconstruct_points_wang(Tdoa, dmax, c); 
        runtime_wang(mi,ri) = toc;
        
         % Error after refinement
        X_aft = align_points(X_c, X_aft, rescale);
        R_aft = X_aft(:, 1:M);
        S_aft = X_aft(:, M+1:end);
        
        err_wang(mi,ri) = (sum(sqrt(sum((R_aft - X_r_c).^2, 1))) + sum(sqrt(sum((S_aft - X_s_c).^2, 1)))) / N;
      
        tic
        [X_bef, X_aft] = reconstruct_points_rand(Tdoa, d, known_delay, 1,1, 'sedumi'); 
        runtime_ours(mi,ri) = toc;
        
%         %% Compute errors
        
        % Error after refinement
        X_aft = align_points(X_c, X_aft, rescale);
        R_aft = X_aft(:, 1:M);
        S_aft = X_aft(:, M+1:end);
        
        err_ours(mi,ri) = (sum(sqrt(sum((R_aft - X_r_c).^2, 1))) + sum(sqrt(sum((S_aft - X_s_c).^2, 1)))) / N;
%         
        
    end
    
    
   
end
% save results
filename = strcat(outDir,'/matlab_', suff, '_M', num2str(M), '.mat');
save(filename,'runtime_ours','runtime_wang', 'err_ours','err_wang')


    