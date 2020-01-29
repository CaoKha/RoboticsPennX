function x_dot = controller(t, state_vec) % x_dot = f(x)
    
    % sprayer offset
    l = 0.1;
    
    % the position and orientation of the robot
    x = state_vec(1);
    y = state_vec(2);
    theta = state_vec(3);

    % your code here
    state_des = desired_path(t);
    x_des = state_des(1);
    y_des = state_des(2);
    x_dot_des = -2*pi*sin(pi*t/5) + pi/2*cos(pi*t/10);
    y_dot_des = pi*cos(pi*t/10) + pi/2*sin(pi*t/10);
%     x_dotdot_des = -2/5*pi^2*cos(pi*t/5) - pi^2/20*sin(pi*t/10);
%     y_dotdot_des = -pi^2/10*sin(pi*t/10) + pi^2/20*cos(pi*t/10);
    k1 = 40;
    k2 = 30;
    Lgh = [cos(theta) -l*sin(theta); sin(theta) l*cos(theta)];
    Lfh = [0; 0];
    input_u = inv(Lgh)*(-Lfh + [x_dot_des + k1*(x_des - x); y_dot_des + k2*(y_des - y)]);
    x_dot = [cos(theta) 0; sin(theta) 0; 0 1]*input_u;
     
end