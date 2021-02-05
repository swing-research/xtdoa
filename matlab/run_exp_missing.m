clear all;
clc

rng(8);

Ms = [8,10,12];
n_runs = 200;
d = 3; % ambient dimension
c = 343; % speed of sound


noises = [1e-3,1e-4,1e-5,1e-6,0]; % noise standard deviation
deletions = linspace(0.01,0.3,10); % percentage of erasures

solver = "sedumi"; % solver fo cvx
known_delay = [0,0]; % sigma, tau
sdr_init = 1; % sdr initialization
rand_restarts = 1;
suff = 'missing';

outDir = strcat('data/',suff);
mkdir(outDir)

rescale = 1;
for mi=1:length(Ms)
    M = Ms(mi); % number of microphones
    K = M; % number of sources
    tot = M*K; % total number of edges
    filename = strcat('data/data_M', num2str(M),'_',num2str(d), 'D.mat');
    data = load(filename);
    data = data.data;
    
    N = M + K;  % number of points
    J_N = eye(N) - 1/N*ones(N);
    err_bef = zeros(n_runs, length(noises), length(deletions)); % error before lm
    err_lm = zeros(n_runs, length(noises),length(deletions)); % final error after refinement

    for si=1:length(noises)
       
        stdev = noises(si);
        rng(8);
        for di=1:length(deletions)
            per = deletions(di); % fraction of missing data
            n = round(per*tot); % number of missing measurements
            
            disp(['Missing data ' num2str(n) ' out of ' num2str(tot)])
            disp([' noise ',num2str(stdev),' (', num2str(M),',' num2str(K),')'])
            for ri = 1:n_runs
                disp(['run ',num2str(ri), ' noise ',num2str(stdev),' (', num2str(M),',' num2str(K),')'])      

               %% Generate points (Forward)
               X = data{ri,4}; % true positions
               X_c = X * J_N; % centered matrices
               X_r_c = X_c(:, 1:M);
               X_s_c = X_c(:, M+1:end);

               T = combine_data(data{ri,1}, data{ri,2}, data{ri,3}, stdev, c, known_delay);

               % create erasures
               rinds = sort(randperm(tot,n)); %pick random positions to delete
               Wm = ones(M,K);
               Wm(rinds) = 0; % mask of missing entries
               Tm = T.*Wm; %matrix with missing entries

               %% Reconstruct points (Inverse)
               [X_bef, X_aft] = reconstruct_points_missing(Tm, d, Wm, known_delay, sdr_init, rand_restarts, solver); % before and after l

               %% Compute errors

               % alignment
               X_bef = align_points(X_c, X_bef, rescale);
               R_bef = X_bef(:, 1:M);
               S_bef = X_bef(:, M+1:end);


               err_bef(ri,si,di) = (sum(sqrt(sum((R_bef - X_r_c).^2, 1))) + sum(sqrt(sum((S_bef - X_s_c).^2, 1)))) / N;

               % Error after refinement
               X_aft = align_points(X_c, X_aft, rescale);       
               R_aft = X_aft(:, 1:M);
               S_aft = X_aft(:, M+1:end);

               err_lm(ri,si,di) = (sum(sqrt(sum((R_aft - X_r_c).^2, 1))) + sum(sqrt(sum((S_aft - X_s_c).^2, 1)))) / N;


            end
        end
    end
    % save results
    filename = strcat(outDir,'/matlab_', suff, '_M', num2str(M), '.mat');
    save(filename,'err_bef','err_lm')
end

    