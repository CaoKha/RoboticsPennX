% number of links to consider
n = 2;

% DH parameters
syms a1 d1 d2 q1 q2 ... 
    Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1 ...
    Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2 real

% the center of mass of each link measured relative to the link fixed frame
% (e.g. c1 = [c1x c1y c1z]' is measured relative to x1y1z1)
c = cell(n,1);
c{1} = [sym('c1x'); sym('c1y'); sym('c1z')];
c{2} = [sym('c2x'); sym('c2y'); sym('c2z')];
assume(vertcat(c{:}), 'real');

% the inertia tensor of each link relative to the inertial frame
I = cell(n,1);
% I{1} = inertia_tensor(1); % please write a function on your own which describe the inertia matrix of the robot you are working with 
% I{2} = inertia_tensor(2); % please write a function on your own which describe the inertia matrix of the robot you are working with
I{1} = [    Ixx1, Ixy1, Ixz1;
            Iyx1, Iyy1, Iyz1;
            Izx1, Izy1, Izz1];
I{2} =  [   Ixx2, Ixy2, Ixz2;
            Iyx2, Iyy2, Iyz2;
            Izx2, Izy2, Izz2];
% mass of each link
syms m1 m2 real

% joint velocities, where qd1 stands for 'q dot 1', the
% first derivative of q1 with respect to time
syms qd1 qd2 real

% acceleration due to gravity (assume g has the correct sign); in other words, if
% gravity were to act in the 'x' direction, the gravity vector would be [g 0 0]
syms g real

% initial conditions for the configuration of Sawyer shown in Figure 1.
% HINT: double(subs(expr, vars, vals)) evaluates a symbolic expression 'expr' by
% substituting each element of 'vals' with its corresponding symbolic variable in 'vars'
q0 = [0 3*pi/2]; % [q10 q20] mm
d0 = [317 192.5]; % [d10 d20] mm
a10 = 81; % in mm

%% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic 
% transformation matrix in terms of the given DH parameters that transforms objects
% in frame i to the inertial frame 0
Ti = cell(n,1);

% The angular velocity Jacobian as an nx1 cell array where each element, Jw{i} is 
% a 3xn symbolic matrix
Jw = cell(n,1);

% The linear velocity Jacobian as an nx1 cell array where each element, Jv{i} is 
% a 3xn symbolic matrix
Jv = cell(n,1);

%% homogeneous transformations
% Transformation from frame 0 to frame 1
Ti{1}=[ cos(q1) -sin(q1)*0  sin(q1)*(-1)    a1*cos(q1)
        sin(q1) cos(q1)*0   -cos(q1)*(-1)   a1*sin(q1)
        0       -1          0               d1
        0       0           0               1];
% Transformation from frame 1 to frame 2    
T12 = [ cos(q2) -sin(q2)*0  sin(q2)*(-1)    0
        sin(q2) cos(q2)*0   -cos(q2)*(-1)   0
        0       -1          0               d2
        0       0           0               1];
% Transformation from frame 0 to frame 2    
Ti{2}= Ti{1}*T12;

%% angular velocity Jacobian
Jw{1} = [[0;0;1] zeros(3,1)];
Jw{2} = [[0;0;1] Ti{1}(1:3,3)];

%% linear velocity Jacobian
% Center of Mass in frame 0
CoM1 = Ti{1}*[c{1};1];
CoM2 = Ti{2}*[c{2};1];
% Jv{n} = [Jv(1) ... Jv(n)], where  Jv(i) = [z(i-1) x (CoM(n) - Origin(i-1))
Jv{1} = [cross([0;0;1],CoM1(1:3)) zeros(3,1)];
Jv{2} = [cross([0;0;1],CoM2(1:3)) cross(Ti{1}(1:3,3),CoM2(1:3) - Ti{1}(1:3,4))];

%% the inertia matrix
% Inertia matrix
D = m1*Jv{1}'*Jv{1} + m2*Jv{2}'*Jv{2} + Jw{1}'*I{1}*Jw{1} + Jw{2}'*I{2}*Jw{2};

% kinetic energy
KE = 1/2*[qd1 qd2]*D*[qd1;qd2];

% potential energy
PE = m1*[0 0 g]*CoM1(1:3) + m2*[0 0 g]*CoM2(1:3);

%% The C matrix 
% Calcualte Chritoffel coefficient matrix
% The grader using the wrong formula in spong.pdf book which is double "dj"
cijk = sym([2 2 2]);
cijk(1,1,1) = 1/2*diff(D(1,1),q1);
cijk(2,1,1) = 1/2*(diff(D(1,1),q1) + diff(D(1,2),q1) - diff(D(2,1),q1));
cijk(1,1,2) = diff(D(2,1),q1) - 1/2*diff(D(1,1),q2);
cijk(2,1,2) = 1/2*(diff(D(2,1),q1) + diff(D(2,2),q1) - diff(D(2,1),q2));
cijk(1,2,1) = 1/2*(diff(D(1,2),q2) + diff(D(1,1),q2) - diff(D(1,2),q1));
cijk(2,2,1) = diff(D(1,2),q2) - 1/2*diff(D(2,2),q1);
cijk(1,2,2) = 1/2*(diff(D(2,2),q2) + diff(D(2,1),q2) - diff(D(1,2),q2)) ;
cijk(2,2,2) = 1/2*diff(D(2,2),q2);

% Correct answer
% cijk(1,1,1) = 1/2*diff(D(1,1),q1);
% cijk(2,1,1) = 1/2*(diff(D(1,1),q2) + diff(D(1,2),q1) - diff(D(2,1),q1));
% cijk(1,1,2) = diff(D(2,1),q1) - 1/2*diff(D(1,1),q2);
% cijk(2,1,2) = 1/2*diff(D(2,2),q1);
% cijk(1,2,1) = 1/2*diff(D(1,1),q2);
% cijk(2,2,1) = diff(D(1,2),q2) - 1/2*diff(D(2,2),q1);
% cijk(1,2,2) = 1/2*(diff(D(2,2),q1) + diff(D(2,1),q2) - diff(D(1,2),q2)) ;
% cijk(2,2,2) = 1/2*diff(D(2,2),q2);

% The grader using Cjk instead of Ckj in the spong.pdf
C = sym([2 2]);
C(1,1) = cijk(1,1,1)*qd1 + cijk(2,1,1)*qd2;
C(1,2) = cijk(1,1,2)*qd1 + cijk(2,1,2)*qd2;
C(2,1) = cijk(1,2,1)*qd1 + cijk(2,2,1)*qd2;
C(2,2) = cijk(1,2,2)*qd1 + cijk(2,2,2)*qd2;

% Correct answer
% C(1,1) = cijk(1,1,1)*qd1 + cijk(2,1,1)*qd2;
% C(1,2) = cijk(1,2,1)*qd1 + cijk(2,2,1)*qd2;
% C(2,1) = cijk(1,1,2)*qd1 + cijk(2,1,2)*qd2;
% C(2,2) = cijk(1,2,2)*qd1 + cijk(2,2,2)*qd2;

% Left side of Lagrange equation of motion 
eom_lhs = D*[qdd1;qdd2] + C*[qd1;qd2] + [diff(PE,q1);diff(PE,q2)];