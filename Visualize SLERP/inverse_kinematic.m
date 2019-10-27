a=13;b=2.5;c=8;d=2.5;e=8;f=2.5;
T01 = compute_dh_matrix(0, pi/2, a, 0);
T12 = compute_dh_matrix(c, 0, -b, 0);
T23 = compute_dh_matrix(0, -pi/2, -d, 0);
T34 = compute_dh_matrix(0, pi/2, e, 0);
T45 = compute_dh_matrix(0, -pi/2, 0, 0);
T56 = compute_dh_matrix(0, 0, f, 0);
T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T36 = T34*T45*T56;




