function [R, S, all_cost] = refine_positions_lm_auglagrange(R0, S0, T, known_delay, inds, S_edm)
% Augmented Lagrangian Algorithm for Equality-constrained LM

M = size(R0, 2);
K = size(S0, 2);
N = M+K;

P = length(inds);

R = R0;
S = S0;
X = [R,S];
Dcurr = edm(X,X);
g = Dcurr(inds) - S_edm; % residual
cost = norm(g);
    
z = zeros(P,1);
mu = 1;

n_iter = 50;
all_cost = zeros(n_iter,1);
thresh = 1e-11;
for i_iter = 1:n_iter
%     mu
    [R,S] = levenberg_marquardt(R,S,T,z,mu,known_delay,inds,S_edm);
    X = [R, S];
    
    % g(x) residual
    Dcurr = edm(X,X);
    g = Dcurr(inds) - S_edm; % residual
    
    z = z + 2*mu*g;
    
    old_cost = cost;
    cost = norm(g);
    all_cost(i_iter) = cost;
    if cost >= 0.25*old_cost 
        mu = 2*mu;
    end
    if cost<thresh
        break
    end
    relch = abs(cost - old_cost) / old_cost;
    if relch < 1e-11
        break
    end
    
    
end


function [R,S] = levenberg_marquardt(R,S,T,z,mu,known_delay,inds,S_edm)

M = size(R, 2);
K = size(S, 2);
d = size(R, 1);
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

I = eye(d *  (M + K));

A = kron(J_K', J_M);
lam = 2;
factor = 0.9;
cost = 1e4;

n_iter = 200;

all_cost = zeros(n_iter,1);
for i_iter = 1:n_iter
    X = [R, S];
    vec_X = X(:);
    
    Delta = sqrt(edm(R, S));
    vec_Delta = Delta(:);

    vec_T = T(:);
    t = A * vec_T;    
    
    % linearization of first term
    Jac = get_jacobian(R, S, known_delay);
    f = A * vec_Delta - t; % residual
  
    % linearization of second term
    JacD = get_jacobian_distance(X, inds);
    Dcurr = edm(X,X);
    g = Dcurr(inds) - S_edm; % residual
    
    old_R = R;
    old_S = S;
    
    vec_X = vec_X - (Jac'*Jac + mu*(JacD'*JacD) + lam*I) \ (Jac'*f + mu*JacD'*(g+ z/(2*mu)));
    
    R = reshape(vec_X(1:d*M), d, M);
    S = reshape(vec_X(d*M+1:end), d, K);
    X = [R, S];
    Dcurr = edm(X,X);
    g = Dcurr(inds) - S_edm; % residual

    old_cost = cost;
    cost = norm(J_M*(sqrt(edm(R, S)) - T)*J_K)^2 + mu*norm(g+z/(2*mu))^2;
    all_cost(i_iter) = cost;
    
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
    end
    
end

function Jac = get_jacobian(R, S, known_delay)
    
M = size(R, 2);
K = size(S, 2);

N = M + K;
d = size(R, 1);

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
                     
function Jac = get_jacobian_distance(X, inds)
    N = size(X,2); % number of points
    d = size(X,1); % dimension of space
    Nd = size(inds,1); % number of known distances
    Jac = zeros(Nd, N*d);
    
    for ni=1:Nd
        curr_ind = inds(ni);
        [indi,indj] = ind2sub([N,N],curr_ind);
        xi = X(:,indi);
        xj = X(:,indj);
        Jac(ni,(indi-1)*d+1:indi*d) = xi - 2*xj;
        Jac(ni,(indj-1)*d+1:indj*d) = -2*xi + xj;
    end
        