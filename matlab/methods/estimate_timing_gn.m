function [sig,eta,X] = estimate_timing_gn(T, dmax)
% Algorithm 1 from Wang et al 2016
% Gauss-Newton algorithm with random restarts
% Estimate the unknown capture and emission times from TDoA matrix
% Input:
%      T: MxK TDoA matrix
%     dmax: maximum distance between sensor and source
%     c: speed of sound
% Output:
%      sig: microphone capture times
%      eta: source emission times
%      X: matrix

[M,K] = size(T); % microphones by sources

% z_T1 = 10^5; % value from paper
z_T1 = 10^10;
z_T2 = 10^-15;
rho_thresh = 10^-9;
gam = 10^9; % weight in objective

rand_restarts = 1000; %max number of initializations
max_iter = 100;% iterations in Gauss-Newton

tij = mean(T,2);

D = T(2:end,2:end).^2 - T(2:end,1).^2- T(1,2:end).^2 + T(1,1).^2;

A = D(:,1:3);
B = D(:,4:end);

% total number of parameters: M + (K-1) + 3*(K-4)
rho0 = zeros(M+K+3*(K-4),1); %output: sig;eta;X
z0 = inf;

for ri=1:rand_restarts
    % initialize
    eta = -dmax + 2*dmax*rand(K-1,1); % emission times
    
    sig = zeros(M,1); % capture times
    sig(1) = (tij(1)-dmax) + 2*dmax*rand(1);
    sig(2:end) = (tij(2:end)-(2*K-1)*dmax/(K)) + 2*dmax*(3*K-2)*rand(M-1,1);
    
    U = get_U(T,eta,sig);
    F = U(:,1:3);
    G = U(:,4:end);
    
    X = (A+F)\(B+G);
    
    rho = [vec(sig);vec(eta);vec(X)]; % variables
    
    for n_iter=1:max_iter
        
        sig = rho(1:M); % M delays
        eta = rho(M+1:M+K-1); % K-1 emissions
        X = reshape(rho(K+M:end),3,K-4);
        old_rho = [vec(sig);vec(eta);vec(X)]; % variables
        
        U = get_U(T,eta,sig);
        F = U(:,1:3);
        G = U(:,4:end);
        
        % current residual
        zA = vec(U);
        zB = vec((A+F)*X - B-G);
        
        Jac = get_jacobian(T,eta,sig,X,gam);
        
        
        % update
        rho = rho - (Jac'*Jac) \ (Jac'*[zA;gam*zB]); % GN
        
        zB_norm = norm((A+F)*X - B- G,'fro')^2;
        
        
        if zB_norm>z_T1 || norm(rho-old_rho)<rho_thresh
            break
        end
        
    end
    
    if zB_norm <z0
%         disp('update rho')
        z0 = zB_norm;
        rho0 = rho;
    end
    
    if z0<z_T2
        break
    end
    
end
sig = rho0(1:M); % M delays
eta = rho0(M+1:M+K-1); % K-1 emissions
X = reshape(rho0(K+M:end),3,K-4);

eta = [0;eta]; % assumption eta(1) = 0

function U = get_U(T,eta,sig)
[M,K] = size(T);
U = zeros(M-1,K-1);
for  mi=2:M
    for ki=2:K
        U(mi-1,ki-1) = 2*sig(mi)*(T(mi,ki)-T(mi,1)) - 2*sig(1)*(T(1,ki)-T(1,1))...
            - 2*eta(ki-1)*(T(mi,ki)-T(1,ki)) + 2*eta(ki-1)*(sig(1)-sig(mi));
    end
end


function Jac = get_jacobian(T, eta, sig, X,lam)

[M,K] = size(T);
P = 3*(K-4);

D = T(2:end,2:end).^2 - T(2:end,1).^2- T(1,2:end).^2 + T(1,1).^2;
U = get_U(T,eta,sig);

A = D(:,1:3);
F = U(:,1:3);
W = A+F;

J11 = zeros(M-1,K-1,M); % sigma derivatives
J12 = zeros(M-1,K-1,K-1); % eta derivatives
J13 = zeros(M-1,K-1,P); % X derivatives


for i=2:M
    for j=2:K
        
        if i==2
            J11(i-1,j-1,1) = -2*(T(1,j)-T(1,1)) + 2*eta(j-1); % sig_1 deriv
        end
        J11(i-1,j-1,i) = 2*(T(i,j)-T(i,1)) - 2*eta(j-1); % sig_i deriv
        
        
%         J12(i-1,j-1,j-1) = -(2*T(i,j)-T(1,j) + 2*sig(i)); % derivation in paper
        J12(i-1,j-1,j-1) = -(2*T(i,j)-T(1,j)) + 2*(sig(1)-sig(i));
        
    end
end

J21 = zeros(M-1,K-4,M); % sigma
J22 = zeros(M-1,K-4,K-1); % eta
J23 = zeros(M-1,K-4,3,K-4); % X

for i=2:M
    for j=5:K
        
        
        for m=1:M
            
            J21(i-1,j-4,m) = J11(i-1,1:3,m)*X(:,j-4) - J11(i-1,j-1,m);
        end
        
        for k=1:K-1
            J22(i-1,j-4,k) = J12(i-1,1:3,k)*X(:,j-4) - J12(i-1,j-1,k);
        end
        
        for k=1:3
            for l=1:K-4
                
                if l==j-4
                    J23(i-1,j-4,k,l) = W(i-1,k);
                end
            end
        end
    end
end

Jac1 = cat(2,reshape(J11,(M-1)*(K-1),M),reshape(J12,(M-1)*(K-1),K-1), reshape(J13,(M-1)*(K-1),P));
Jac2 = cat(2,lam*reshape(J21,(M-1)*(K-4),M),lam*reshape(J22,(M-1)*(K-4),K-1),lam*reshape(J23,(M-1)*(K-4),3*(K-4)));
Jac = cat(1,Jac1,Jac2);

