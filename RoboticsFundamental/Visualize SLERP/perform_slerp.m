%% Helper script to visualise SLERP.
% Remember to fill in your implementation in the quat_slerp function before
% you run this script!

close all
pause on;  % Set this to on if you want to watch the animation
GraphingTimeDelay = 0.05; % The length of time that Matlab should pause between positions when graphing, if at all, in seconds.

% Generate the starting and final orientation
% Start with something simple like this, but feel free to play around with
% different values for q0 and q1!
% Remember to test the case with theta > 180, to ensure the shortest path
% is always chosen
q0 = [1, 0, 0, 0];
theta = 7*pi/4;
q1 = [cos(theta/2), sin(theta/2), 0, 0];
% Once that begins to work, you can try something a little more complex
% like below
% q1 = [cos(theta/2), (1/sqrt(3))*sin(theta/2), (1/sqrt(3))*sin(theta/2), (1/sqrt(3))*sin(theta/2)];

% The 3 axis
x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Setup the plot
figure
scale_f = 2;
axis(scale_f*[-1 1 -1 1 -1 1])
axis vis3d
grid on
view(125,25)
xlabel('X')
ylabel('Y')
zlabel('Z')

% Get the axis of starting frame
x0 = quat_multiply(quat_multiply(q0, [0 x]), quat_conjugate(q0));
y0 = quat_multiply(quat_multiply(q0, [0 y]), quat_conjugate(q0));
z0 = quat_multiply(quat_multiply(q0, [0 z]), quat_conjugate(q0));

o = [0 0 0];
x0 = x0(2:4);
y0 = y0(2:4);
z0 = z0(2:4);

% Get the axis of final frame
x1 = quat_multiply(quat_multiply(q1, [0 x]), quat_conjugate(q1));
y1 = quat_multiply(quat_multiply(q1, [0 y]), quat_conjugate(q1));
z1 = quat_multiply(quat_multiply(q1, [0 z]), quat_conjugate(q1));

x1 = x1(2:4);
y1 = y1(2:4);
z1 = z1(2:4);

% Plot axis
hold on

plot_line(o, x0, 'r');
plot_line(o, y0, 'g');
plot_line(o, z0, 'b');

plot_line(o, x1, 'r--');
plot_line(o, y1, 'g--');
plot_line(o, z1, 'b--');

% Animate the vector
q_anim = quat_slerp( q0, q1, 50 );

hx = plot_line(o, x0, 'r--');
hy = plot_line(o, y0, 'g--');
hz = plot_line(o, z0, 'b--');

pause(GraphingTimeDelay);
for i = 1:size(q_anim, 1)
   
    % Plot vector
    v_animX = quat_multiply(quat_multiply(q_anim(i, :), [0 x]), quat_conjugate(q_anim(i, :)));
    v_animY = quat_multiply(quat_multiply(q_anim(i, :), [0 y]), quat_conjugate(q_anim(i, :)));
    v_animZ = quat_multiply(quat_multiply(q_anim(i, :), [0 z]), quat_conjugate(q_anim(i, :)));
    v_animX = v_animX(2:4);
    v_animY = v_animY(2:4);
    v_animZ = v_animZ(2:4);
    
    set(hx,'xdata',[0 v_animX(1)], 'ydata', [0 v_animX(2)], 'zdata', [0 v_animX(3)]);
    set(hy,'xdata',[0 v_animY(1)], 'ydata', [0 v_animY(2)], 'zdata', [0 v_animY(3)]);
    set(hz,'xdata',[0 v_animZ(1)], 'ydata', [0 v_animZ(2)], 'zdata', [0 v_animZ(3)]);
    drawnow
    
    pause(GraphingTimeDelay);
    
end