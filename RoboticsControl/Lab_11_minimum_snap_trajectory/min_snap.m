function coefficients = min_snap(T,boundary)
    % boundary is a column vector in the form
    % [r0; v0; a0; j0; rT; vT; aT; jT]
    % representing the derivatives of the trajectory
    % at time 0 and time T
    
    syms D(t) 
    D(t) = derivative_matrix(4);
    
    % your implementation...
    D_init = subs(D(t),t,0);
    D_final = subs(D(t),t,T);
    big_D = [D_init ; D_final];
    coefficients = inv(big_D)*boundary;

end