function [X2p, R, t, s] = align_points(X1,X2,rescale)
m1 = mean(X1,2);
m2 = mean(X2,2);
[U,~,V] = svd((X1-m1)*(X2-m2)');
R = U*V';

if(rescale)
   s = norm(X1-m1) / norm(X2-m2);
else
   s = 1.0; 
end

t = m1 - s*R*m2;
X2p = s*R*X2 + t;
end

