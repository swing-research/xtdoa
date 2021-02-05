clear all;
% Script to generate and save dataset

rng(123);

d = 3;          % ambient dimension
runs = 200;


suff = 'known_offset';

for M=6:15  % number of microphones
    K = M; % number of sources
    
    data = cell(runs,4);
    for ri=1:runs
        [T,sig,tau,X] = generate_points(d,M,K);
        data{ri,1} = T;
        data{ri,2} = sig;
        data{ri,3} = tau;
        data{ri,4} = X;
    end
    filename = strcat('data/data_M', num2str(M),'_',num2str(d), 'D.mat');
    save(filename,'data')
end
    
