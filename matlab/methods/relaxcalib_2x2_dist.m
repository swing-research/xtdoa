function X_est = relaxcalib_2x2_dist(T, d, known_delay, solver, inds, S_edm)


M = size(T, 1);
K = size(T, 2);
N = M + K;

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

S_row = [eye(M) zeros(M, K)];
S_col = [zeros(M, K); eye(K)];


cvx_begin sdp quiet
    cvx_precision low
    if solver=="sdpt3"
        cvx_solver sdpt3
    elseif solver=="sedumi"
        cvx_solver sedumi
    end

    variable G(N, N) symmetric;
    variable B(M, K);

    DofG = diag(G)*ones(N, 1)' - 2*G + ones(N, 1)*diag(G)';
    LofG = S_row*DofG*S_col;

    G >= 0;
    G*ones(N, 1) == 0;
    DofG(inds)==S_edm; % known distances

    L = cell(M, K);
    for m = 1:M
        for k = 1:K
            L{m, k} = [LofG(m, k) B(m, k); B(m, k) 1];
            L{m, k} >= 0;
        end
    end

    B(:) >= 0;

    minimize square_pos(norm(J_M*(B - T)*J_K, 'fro')) %+ square_pos(norm(DofG(inds)-S_edm))
cvx_end


[~, S, V] = svd(G);
X_est = sqrt(S(1:d, :))*V';

end