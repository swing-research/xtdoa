function [X_est, A] = relaxcalib_2x2_missing(T, d, W, known_delay)

M = size(T, 1);
K = size(T, 2);
N = M + K;

J_M = eye(M) - 1/M*ones(M);
J_K = eye(K) - 1/K*ones(K);
J_N = eye(N) - 1/N*ones(N);

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


num = sum(W==0,'all'); % number of missing measurements
[indm,indk] = find(W==0); %indices of missing measurements
inds = find(W==0);

e = @(k,n) [zeros(k-1,1);1;zeros(n-k,1)]; % standard n-dim basis
% E = @(alpha) sum(alpha(i)*e(indm(i),M)*e(indk(i),K)');

S_row = [eye(M) zeros(M, K)];
S_col = [zeros(M, K); eye(K)];
const= 1e-20;

cvx_begin sdp quiet
    cvx_precision low
%     cvx_solver sdpt3
    cvx_solver sedumi

    variable G(N, N) symmetric;
    variable B(M, K);
    variable A(num,1);
%     expression E(M,K);
     E = cvx(zeros(M,K));

    DofG = S_row*(diag(G)*ones(N, 1)' - 2*G + ones(N, 1)*diag(G)')*S_col;

    G >= 0;
    G*ones(N, 1) == 0;

    L = cell(M, K);
    for m = 1:M
        for k = 1:K
            L{m, k} = [DofG(m, k) B(m, k); B(m, k) 1];
            L{m, k} >= 0;
        end
    end

    B(:) >= 0;
%     E = zeros(M,K);
    E(inds) = A;
%     for ai=1:num
%         E = E +  A(ai)*e(indm(ai),M)*e(indk(ai),K)';
%     end
        minimize square_pos(norm(J_M*(B + E - T)*J_K, 'fro'))
%         % handling missing part with alpha
cvx_end
     
[~, S, V] = svd(G);
X_est = sqrt(S(1:d, :))*V';

% if ~any(isnan(G)) & ~any(isinf(G))
% [~, S, V] = svd(G);
% X_est = sqrt(S(1:d, :))*V';
% else
%     X_est = zeros(d,N);
% end

% errBT = norm(B-T,'fro');
% display(['error B-T ' num2str(errBT)])
% figure;
% imagesc(B)
% colorbar()
% title('B')

% figure;
% imagesc(T)
% colorbar()
% title('T')


end
