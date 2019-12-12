clear; clc;

%% homogeneous transforms

n = 7; % DOF

% DH parameters
q = sym('q', [n 1], 'real'); % generalized coordinates (joint angles)
d = sym('d', [n 1], 'real'); % link offsets
syms a1 real

% initial conditions for the configuration of Sawyer shown in Figure 1.
% you can use these values to sense check your work
% HINT: vpa(subs(expr, vars, vals)) evaluates a symbolic expression 'expr' by
% substituting each element of 'vals' with its corresponding symbolic variable in 'vars'
q0 = [0 3*pi/2 0 pi 0 pi 3*pi/2];
d0 = [317 192.5 400 168.5 400 136.3 133.75];
a10 = 81;

% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic transform matrix
% provide your answer in terms of the given symbolic variables
% NOTE: for symbolic arrays: q(1) = q1, q(2) = q2, etc.
Ti = cell(n,1);

% Initialize array of alpha (rotation by x axis)
alpha = sym([n 1]);
alpha(1) = -pi/2;
alpha(2) = -pi/2;
alpha(3) = -pi/2;
alpha(4) = -pi/2;
alpha(5) = -pi/2;
alpha(6) = -pi/2;
alpha(7) = 0;

%Initialize array of distance a (zi-1 to zi via xi)
a = sym([n 1]);
a(1) = a1;
a(2:n) = 0;

% Ti = your homogeneous transformations solution
for i = 1:n
    % input symbolic homogeneous transformation matrix i
    % Calculate T[(i-1)->i]
    Ti{i} = [   cos(q(i)), -sin(q(i))*cos(alpha(i)), sin(q(i))*sin(alpha(i)), a(i)*cos(q(i));
                sin(q(i)), cos(q(i))*cos(alpha(i)), -cos(q(i))*sin(alpha(i)), a(i)*sin(q(i));
                0,         sin(alpha(i)),           cos(alpha(i)),            d(i);
                0,          0,                      0,                        1];
    % Calculate T[0->i]
    if i ~= 1
        Ti{i} = Ti{i-1}*Ti{i};
    end
end
%% angular velocity jacobian (Jw)

% Initialize angular velocity jacobian as an nx1 cell array where each element is
% an 3xn symbolic matrix
Jw = arrayfun(@(x) sym(['Jw' num2str(x)], [3,n], 'real'), 1:n, 'UniformOutput', 0)';
% Jw = your angular velocity jacobian solution

% Initialize array z
z = cell(n,1);
for i = 1:n
    z{i} = Ti{i}(1:3,3);
end

for i = 1:n
    Jw{i} = sym([3 n]);
end

% Calculate Jw
for i = 1:n
    for k = 1:i
        if k == 1
            Jw{i}(1:3,k) = [0;0;1];
        else
            Jw{i}(1:3,k) = z{k-1};
        end
    end
    for k = i+1:n
        if (i+1)<=n
            Jw{i}(1:3,k) = [0;0;0];
        end
    end
end
%% linear velocity jacobian (Jv)

% the center of mass of each link measured relative to the link fixed frame
% like Ti and Jw, c is an nx1 cell array where each element is a symoblic vector/matrix
% for example: c{3} = [c3x c3y c3z]' is the center of mass of link 3 measured relative to frame 3
c = arrayfun(@(x) [sym(['c' num2str(x) 'x'], 'real'), sym(['c' num2str(x) 'y'], 'real'), ...
    sym(['c' num2str(x) 'z'], 'real')]', 1:n, 'UniformOutput', 0)';

% as with the angular velocity jacobian, the linear velocity jacobian is comprised of n 3xn
% symbolic matrices stored in a cell array. Jv{i} is the 3xn angular velocity jacobian of link i
Jv = cell(n,1);

%  Jv = your linear velocity jacobian solution
for i = 1:n
    Jv{i} = sym([3 n]);
end

CoM = cell(n,1);

for i = 1:n
    CoM{i} = Ti{i}*[c{i};1];
end
% Calculate Jv
for i = 1:n
    for k = 1:i
        if k == 1
            Jv{i}(1:3,k) = cross([0;0;1],CoM{i}(1:3));
        else
            Jv{i}(1:3,k) = cross(z{k-1},CoM{i}(1:3) - Ti{k-1}(1:3,4));
        end
    end
    for k = i+1:n
        if (i+1)<=n
            Jv{i}(1:3,k) = [0;0;0];
        end
    end
end
%% potential energy

m = sym('m', [n 1], 'real'); % mass of each link
syms g real % gravity
% PE = your potential energy solution
PE = 0;
for i = 1:n
    PE = PE + m(i)*[0 0 g]*CoM{i}(1:3);% total potential energy of Sawyer
end

%% inertial matrix and kinetic energy

qd = sym('qd', [n 1], 'real'); % "q dot" - the first derivative of the q's in time (joint velocities)

% inertia tensor for each link relative to the inertial frame stored in an nx1 cell array
I = arrayfun(@(x) inertia_tensor(x), 1:n, 'UniformOutput', 0)';

% D = your inertia matrix solution
D = 0;
for i = 1:n
    D = D + m(i)*Jv{i}'*Jv{i} + Jw{i}'*I{i}*Jw{i}; 
end
% KE = your kinetic energy solution
KE = 1/2*qd'*D*qd;
%% equations of motion

qdd = sym('qdd', [n 1], 'real'); % "q double dot" - the second derivative of the q's in time (joint accelerations)
% Compute C
% cijk = sym([n n n]);
% for i = 1:n
%     for j = 1:n
%         for k = 1:n
%             cijk(i,j,k) = 1/2*(diff(D(k,j),q(j))+diff(D(k,i),q(j))-diff(D(i,j),q(k)));
%         end 
%     end
% end

C = sym(zeros(n, n));
for j = 1:n
    for k = 1:n
        for i = 1:n
            C(j,k) = C(j,k) + 1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qd(i);
        end
    end
end
% Compute gk
gk = sym(zeros(n,1));
for i = 1:n
    gk(i,1) = diff(PE,q(i));
end
% eom_lhs = your solution
eom_lhs = D*qdd + C*qd + gk;
%%%%%%% THIS IS THE END OF YOUR INPUT/EDITS %%%%%%%%

%% helper functions (don't use/edit - only used to help initialize I)
function tensor = inertia_tensor(num)

n = num2str(num);

tensor = [sym(['Ixx' n]) sym(['Ixy' n]) sym(['Ixz' n]);
          sym(['Iyx' n]) sym(['Iyy' n]) sym(['Iyz' n]);
          sym(['Izx' n]) sym(['Izy' n]) sym(['Izz' n])];

assume(tensor, 'real');

end