function [R, S, cost] = refine_positions_lm(R0, S0, T, known_delay)
% Levenberg-Marquardt Algorithm for timing-invariant localization from
% times-of-arrival
% Input:
%      T: MxK ToA or TDoA matrix
%      R0: dxM initial receiver positions
%      S0: dxK initial source positions
%     known_delay: list indicating if sensors or sources are synchronized
%     [0,0] if fully unsynchronized
% Output:
%      R: dxM receiver position reconstruction
%      S: dxK source point reconstruction
%      cost: final cost


M = size(R0, 2);
K = size(S0, 2);
J_M = eye(M) - (1/M)*ones(M);
J_K = eye(K) - (1/K)*ones(K);

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

d = size(R0, 1);
I = eye(d *  (M + K));

A = kron(J_K', J_M);
R = R0;
S = S0;

lam = 2;
lam_min = 0.1;
factor = 0.9;
cost = 1e4;

n_iter = 1000;

for i_iter = 1:n_iter
    X = [R, S];
    vec_X = X(:);
    
    Delta = sqrt(edm(R, S));
    vec_Delta = Delta(:);

    vec_T = T(:);
    t = A * vec_T;    
    
    Jac = get_jacobian(R, S, known_delay);
    
    old_R = R;
    old_S = S;
    
    f = A * vec_Delta - t;
    vec_X = vec_X - (Jac' * Jac + lam * I) \ (Jac' * f);
    
    R = reshape(vec_X(1:d*M), d, M);
    S = reshape(vec_X(d*M+1:end), d, K);

    old_cost = cost;
    cost = norm(J_M*(sqrt(edm(R, S)) - T)*J_K);
    
    
    if cost <= old_cost
        lam = factor * lam;
        relch = abs(cost - old_cost) / old_cost;
        if relch < 1e-11
%             fprintf('Relative tolerance reached\n');
            break
        end
    else
        lam = 2 * lam;
        cost = old_cost;
        R = old_R;
        S = old_S;
    end
end


function Jac = get_jacobian(R, S, known_delay)
    
M = size(R, 2);
K = size(S, 2);

N = M + K;
d = size(R, 1);


J_M = eye(M) - (1/M)*ones(M);
J_K = eye(K) - (1/K)*ones(K);
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

for m = 1:M
    for k = 1:K            

        mictosrc = R(:, m) - S(:, k);
        JacD_R(m, k, :, m) = mictosrc / norm(mictosrc);
        JacD_S(m, k, :, k) = -mictosrc / norm(mictosrc);
    end
end

JacD = cat(4, JacD_R, JacD_S);
JacD_mtx = reshape(JacD, M*K, d*N); 
A = kron(J_K', J_M);
Jac = A * JacD_mtx; 