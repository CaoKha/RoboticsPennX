function [ q_int ] = quat_slerp( q0, q1, steps )
%QUAT_SLERP Perform SLERP between two quaternions and return the intermediate quaternions
%   Usage: [ q_int ] = quat_slerp( q0, q1, steps )
%   Inputs:
%       q0 is the quaternion representing the starting orientation, 1x4 matrix
%       q1 is the quaternion representing the final orientation, 1x4 matrix
%       steps is the number of intermediate quaternions required to be returned, integer value
%       The first step is q0, and the last step is q1
%   Output:
%       q_int contains q0, steps-2 intermediate quaternions, q1
%       q_int is a (steps x 4) matrix

    %% Your code goes here
    q_int = zeros(steps, 4);
    q_int(1,:) = q0;
    q_int(steps,:) = q1;
    cosHalfTheta = dot(q0,q1);     
    if (cosHalfTheta < 0) 
        q1 = -q1;
        cosHalfTheta = -cosHalfTheta;
    end
    halfTheta = acos(cosHalfTheta);
    k = 1;
    for i = 0:(1/(steps-1)):1
        q_int(k,:)= (sin((1-i)*halfTheta))/sin(halfTheta)*q0 + sin(i*halfTheta)/sin(halfTheta)*q1;
        k = k+ 1;
        
    end 
end