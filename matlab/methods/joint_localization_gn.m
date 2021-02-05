function [R,S,beta0] = joint_localization_gn(T,sig,eta)
% Algorithm 2 from Wang et al 2016
% Gauss-Newton algorithm with random restarts
% Estimate the source and sensor positions from TDoA matrix

[M,K] = size(T); % microphones by sources

assert(length(eta)==K)

% v_T1 = 10^5; % value from paper, doesn't work
v_T1 = 10^10;
v_T2 = 10^-15;
beta_thresh = 10^-9;


rand_restarts = 50; %max number of initializations
max_iter = 100;% iterations in Gauss-Newton


D = T(2:end,2:end).^2 - T(2:end,1).^2- T(1,2:end).^2 + T(1,1).^2;
U = get_U(T,eta,sig);

D_hat = D+U;

Gij = zeros(M,K);
for m=1:M
    for k=1:K
        Gij(m,k) = (T(m,k)^2 + eta(k)^2 + sig(m)^2 - 2*(T(m,k)*eta(k)-T(m,k)*sig(m)+eta(k)*sig(m)));
    end
end
g_dash = Gij(2:end,2:end)- Gij(1,2:end) - D_hat;
g_dash = g_dash(:,1);
g_ddash = Gij(2:end,2:end)- Gij(2:end,1) - D_hat;
g_ddash = g_ddash(1,:)';

RtS = -1/2 *D_hat;
[D_l,Sig3,D_r] = svd(RtS);
D_l = D_l(:,1:3);
Sig3 = Sig3(1:3,1:3);
D_r = D_r(:,1:3);
    
P = 9+1;
beta0 = zeros(P,1); %variables: sig;eta;X

v0 = inf;

for ri=1:rand_restarts
    % initialize
    beta = randn(P,1);
   
   costall = zeros(max_iter,1);
    for n_iter=1:max_iter

        C = reshape(beta(1:9),3,3);
        sx1 = beta(10:end);
        old_beta = [vec(C);sx1];
        
        R_bar = (D_l*C)';
        S_bar = C\(Sig3*D_r');     
    
        % current residual
        vA = diag(R_bar'*R_bar) - 2*R_bar(1,:)'*sx1 - g_dash;
        vB = diag(S_bar'*S_bar) + 2*S_bar(1,:)'*sx1 - g_ddash;
         
        Jac = get_jacobian(C, sx1, RtS);

        % update
        beta = beta - (Jac' * Jac) \ (Jac' * [vA;vB]); % GN
        
        cost = norm(vA)^2 + norm(vB)^2;
    
        costall(n_iter) = norm(vA)^2 + norm(vB)^2;
        
        if cost>v_T1 || norm(beta-old_beta)<beta_thresh
%             disp('divergence or small change')
%             cost
%             norm(beta-old_beta)
            
            break
        end

    end
    if cost<v0 
%         disp('update beta')
        v0 =cost;
        beta0 = beta;
    end
    
    if v0<v_T2
%         v0
%         disp('global convergence threshold')
        break
    end
   
   
end

C = reshape(beta0(1:9),3,3); %
sx1 = beta0(10); %

R_bar = (D_l*C)';
S_bar = C\(Sig3*D_r');

s1 = [sx1;0;0];
R = [zeros(3,1),R_bar]; % first sensor at origin (0,0,0)
S = [s1,S_bar+s1]; % first receiver at (sx1,0,0)

function U = get_U(T,eta,sig)
[M,K] = size(T);
U = zeros(M-1,K-1);
for m=2:M
    for k=2:K
        U(m-1,k-1) = 2*sig(m)*(T(m,k)-T(m,1)) - 2*sig(1)*(T(1,k)-T(1,1))...
            - 2*eta(k)*(T(m,k)-T(1,k)) + 2*eta(k)*(sig(1)-sig(m)); 
    end
end


function Jac = get_jacobian(C, sx1,RtS)

[D_l,Sig3,D_r] = svd(RtS);
D_l = D_l(:,1:3);
Sig3 = Sig3(1:3,1:3);
D_r = D_r(:,1:3);
    
R_bar = (D_l*C)';
S_bar = C\(Sig3*D_r');

M = size(R_bar,2) + 1;
K = size(S_bar,2) + 1;
        
J11 = zeros(M-1,3,3);
J12 = zeros(M-1,1);

J21 = zeros(K-1,3,3);
J22 = zeros(K-1,1);
B = Sig3*D_r';
% size(D_l)
for k=1:3
    for l=1:3
        

        for i=1:M-1
            J11(i,k,l) = 2*R_bar(l,i)*D_l(i,k);
            if l==1
                J11(i,k,l) = 2*R_bar(l,i)*D_l(i,k) - 2*D_l(i,k)*sx1;
            end
            J12(i,1) = -2*R_bar(1,i);
        end

        Okl = zeros(3,3);
        Okl(k,l) = 1;
        CC = -C\Okl/C;
        S_deriv = CC*B; % derivative of S wrt C_kl
        for j=1:K-1
            
            J21(j,k,l) = 2*S_bar(:,j)'*S_deriv(:,j) + 2*sx1*S_deriv(1,j);
            J22(j,1) = 2*S_bar(1,j);
            
        end
          
    end
    
end

Jac1 = cat(2,reshape(J11,M-1,9),J12);
Jac2 = cat(2,reshape(J21,K-1,9),J22);
Jac = cat(1,Jac1,Jac2);
