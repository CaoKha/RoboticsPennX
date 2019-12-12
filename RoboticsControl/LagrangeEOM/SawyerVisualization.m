q = [0 -pi/2 0 -pi 0 -pi 0];
sawyer_coordinates(q);
T = cell([7 1]);
for i=1:7
    T{i} = [eye(3) [0;0;.2]; zeros(1,3) 1];
end
sawyer_transforms(T);
q = [pi/5 -pi/2 -pi/4 3*pi/4 3*pi/4 -pi/4 pi/8];
qdot = [.015 -.015 .02 -.02 -.05 0 .05];
for i=1:100
    sawyer_coordinates(q+i*qdot);
end
