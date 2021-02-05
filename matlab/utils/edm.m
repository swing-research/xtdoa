function D = edm(R, S)

    M = size(R, 2);
    K = size(S, 2);

    norm2_R = sum(R.^2, 1)';
    norm2_S = sum(S.^2, 1)';

    D = norm2_R*ones(K, 1)' - 2*R'*S + ones(M, 1)*norm2_S';
