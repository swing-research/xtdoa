function [R, S, A, cost] = refine_positions_lm_missing(R0, S0, A0, T, W, known_delay)

M = size(R0, 2);
K = size(S0, 2);

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

num = sum(W==0, 'all'); % number of missing measurements
[indm,indk] = find(W==0); %indices of missing measurements

e = @(k,n) [zeros(k-1,1);1;zeros(n-k,1)]; % standard n-dim basis

d = size(R0, 1);
I = eye(d *  (M + K)+num);

R = R0;
S = S0;
A = A0;

E = zeros(M,K);
for ai=1:num
     E = E +  A(ai)*e(indm(ai),M)*e(indk(ai),K)';
end
    
lam = 2;
lam_min = 0.1;
factor = 0.9;
cost = 1e4;

n_iter = 100;

for i_iter = 1:n_iter

    vec_X = [R(:); S(:); A];
    
    old_R = R;
    old_S = S;
    old_A = A;
    old_E = E;
    
    Jac = get_jacobian(R, S, W, known_delay);
    
    f = fres(R,S,E,T);
    vec_X = vec_X - (Jac' * Jac + lam * I) \ (Jac' * f);
    
    R = reshape(vec_X(1:d*M), d, M);
    S = reshape(vec_X(d*M+1:d*M+d*K), d, K);
    A = vec_X(d*(M+K)+1:end);
    
    E = zeros(M,K);
    for ai=1:num
         E = E +  A(ai)*e(indm(ai),M)*e(indk(ai),K)';
    end

    old_cost = cost;
    cost = norm(J_M*(sqrt(edm(R, S)) + E - T)*J_K);
    
    
    if cost <= old_cost
        lam = factor * lam;
        relch = abs(cost - old_cost) / old_cost;
        if relch < 1e-11
            break
        end
    else
        lam = 2 * lam;
        cost = old_cost;
        R = old_R;
        S = old_S;
        A = old_A;
        E = old_E;
    end
    
end

function fcurr = fres(R,S,E,T)
    M = size(R, 2);
    K = size(S, 2);

    J_M = eye(M) - (1/M)*ones(M);
    J_K = eye(K) - (1/K)*ones(K);

    fcurr = J_M*(sqrt(edm(R, S)) + E - T)*J_K;
    fcurr = fcurr(:);

function Jac = get_jacobian(R, S, W, known_delay)
    
M = size(R, 2);
K = size(S, 2);
d = size(R, 1);

num = sum(W==0,'all'); % number of missing measurements
inds = find(W==0); %indices of missing measurements
[indm,indk] = find(W==0); %indices of missing measurements
% N = d*M + d*K + M*K; % total number of variables

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

JacD_R = zeros(M, K, d, M);
JacD_S = zeros(M, K, d, K);
JacD_E = zeros(M, K, num);
for m = 1:M
    for k = 1:K            
       
        mictosrc = R(:, m) - S(:, k);
        
        JacD_R(m, k, :, m) = mictosrc / norm(mictosrc);
        JacD_S(m, k, :, k) = -mictosrc / norm(mictosrc);
        linind = sub2ind([M,K],m,k);
        if ismember(linind,inds)
            pos = find(inds==linind);
            JacD_E(m, k, pos) = 1;
        end
       
    end
end

JacD = cat(4, JacD_R, JacD_S);
JacD_mtx = reshape(JacD, M*K, d*(M+K));
JacE_mtx = reshape(JacD_E, M*K, num); 
JacD_mtx = cat(2, JacD_mtx, JacE_mtx);

A = kron(J_K', J_M);

Jac = A * JacD_mtx; 