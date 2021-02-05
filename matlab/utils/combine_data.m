function T = combine_data(T, sig, tau, stdev, c, known_delay)

[M,K] = size(T);
% sig and tau are swapped compared to paper!
if ~known_delay(1)
    T = T + c*(ones(M, 1)*sig');
end

if ~known_delay(2)
    T = T + c*(tau*ones(K, 1)');
end

% Add noise
T = T + c*stdev*randn(M, K);