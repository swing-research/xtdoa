function [T,sig,tau,X] = generate_points(d, M, K)

if d == 2 
    roomdim = [10, 7];
else
    roomdim = [10, 10, 3];
end
    
X_r = diag(roomdim) * rand(d, M);
X_s = diag(roomdim) * rand(d, K);
X = [X_r, X_s];


% Random time offsets and internal delays
tau = 2*rand(M, 1) - 1;
sig = 2*rand(K, 1) - 1;


% Matrix of times of arrival and the ground truth (squared) EDM
T = zeros(M, K);
for i = 1:d
    T = T + bsxfun(@minus, X_r(i, :)', X_s(i, :)).^2;
end

T = sqrt(T); 
