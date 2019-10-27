obstacle = boxFV(10, 20, 10, 20);
obstacle = appendFV (obstacle, boxFV(-20, 0, -20, -10));
obstacle = appendFV (obstacle, transformFV(boxFV(-10, 10, -10, 10), 30, [-20 20]));

x1 = [0; 0; 0; 0; 0; 0];
x2 = [30; 30; 30; 30; 30; 30];

out = LocalPlannerSixLink (x1, x2, obstacle)