[X,Y] = meshgrid(-2:0.75:2);
R = sqrt(X.^2 + Y.^2)+ eps;
V = sin(R)./(R);

figure
surf(X,Y,V)
xlim([-4 4])
ylim([-4 4])
title('Original Sampling')