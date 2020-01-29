syms p a B real % p = rho, a = alpha, b = beta 
syms k_p k_a k_b real
% x_dot = f(x), Linearize f(x) at [p a B] = [ 0 0 0 ], f(x) = f(0) + df/dx(x=0)
A = [-k_p 0 0; 0 k_p-k_a -k_b; 0 -k_p 0];% linearized system
% we want all the eigenvalues of A to be real negative
% use eig(A) to find these eigenvalues
% we found that kp > 0 then pick kp = 5, pick randomly ka (ka = 10) then -1.25 < kb < 0 
K = [5 10 -1.24];% stable gains [k_p k_a k_b];
% let's test it out
eig(subs(A,[k_p,k_a,k_b],K))