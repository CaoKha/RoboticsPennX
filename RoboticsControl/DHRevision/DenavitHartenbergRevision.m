syms a1 d1 d2 q1 q2

% provide your answer in terms of the symbolic variables above
R01 = [ cos(q1) -sin(q1)*0  sin(q1)*(-1)
        sin(q1) cos(q1)*0   -cos(q1)*(-1)
        0       -1          0];

R12 = [ cos(q2) -sin(q2)*0  sin(q2)*(-1)
        sin(q2) cos(q2)*0   -cos(q2)*(-1)
        0       -1          0];

R02 = R01*R12;

T01 = [ cos(q1) -sin(q1)*0  sin(q1)*(-1)    a1*cos(q1)
        sin(q1) cos(q1)*0   -cos(q1)*(-1)   a1*sin(q1)
        0       -1          0               d1
        0       0           0               1];

T12 = [ cos(q2) -sin(q2)*0  sin(q2)*(-1)    0
        sin(q2) cos(q2)*0   -cos(q2)*(-1)   0
        0       -1          0               d2
        0       0           0               1];

T02 = T01*T12;

T20 = [ R02'        -R02'*T02(1:3,4)
        0   0   0   1];