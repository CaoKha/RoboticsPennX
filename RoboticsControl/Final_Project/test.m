c = min_snap(5,[0; 0; 0; 0; 3; 0; 0; 0]);
t = 5;
plot_traj(c,t);

boundary = [0; 0; 0; 0; 1; 0; 0; 0];
ps = [2; 3];
ts = [1; 3; 5];
waypoints(boundary, ps, ts)