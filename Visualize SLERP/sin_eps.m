function A = sin_eps (X)
% Compute the cosine for each element of X in radians.
% It returns zero where output is less than epsilon.
%
% For example, sum(cos_eps(pi*(0:0.5:10000))==0)
% It returns
% 	ans = 10000
%
% See also cos, cosd
A = sin(X);
idx_less_eps = abs(A)<eps(X);
A(idx_less_eps) = 0;
return