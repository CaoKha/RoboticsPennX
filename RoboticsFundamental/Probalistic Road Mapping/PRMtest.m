fv = SixLinkRobot ([-120 120 -120 120 -120 60]);

fv2 = SixLinkRobot ([0 0 0 0 0 180]);

obstacle = boxFV(10, 20, 10, 20);
obstacle = appendFV (obstacle, boxFV(-20, 0, -20, -10));
obstacle = appendFV (obstacle, transformFV(boxFV(-10, 10, -10, 10), 30, [-20 20]));

nsamples = 20;
neighbors = 5;
Dist = @DistSixLink;
LocalPlanner = @(x,y)(LocalPlannerSixLink(x,y,obstacle));

roadmap = PRM (@()(RandomSampleSixLink(obstacle)), @DistSixLink, @(x,y)(LocalPlannerSixLink(x,y,obstacle)), nsamples, neighbors);