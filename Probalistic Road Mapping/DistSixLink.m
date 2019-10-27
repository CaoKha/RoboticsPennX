function out = DistSixLink (x1, x2)
% Compute the distance between two sets of six link coordinates
% Note we assume all angular coordinates are between 0 and 360
% Use sum of the absolute value of angle difference (L1 norm) as the
% distance.
% Note this is angle difference.
% i.e. DistSixLink([0, 0, 0, 0, 0, 0], [360, 360, 360, 360, 360, 360]) = 0
e = abs(bsxfun(@minus, x2, x1));
out = sum(min(e, 360-e));
