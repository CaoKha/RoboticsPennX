function A = compute_dh_matrix(r, alpha, d, theta)
%% It might be helpful to copy over this function from Week 5 Lab
    A = zeros(4);
    A(1,1) = cos(theta);
    A(1,2) = -sin(theta)*cos(alpha);
    A(1,3) = sin(theta)*sin(alpha);
    A(1,4) = r*cos(theta);
    A(2,1) = sin(theta);
    A(2,2) = cos(theta)*cos(alpha);
    A(2,3) = -cos(theta)*sin(alpha);
    A(2,4) = r*sin(theta);
    A(3,1) = 0;
    A(3,2) = sin(alpha);
    A(3,3) = cos(alpha);
    A(3,4) = d;
    A(4,1) = 0;
    A(4,2) = 0;
    A(4,3) = 0;
    A(4,4) = 1;
end

