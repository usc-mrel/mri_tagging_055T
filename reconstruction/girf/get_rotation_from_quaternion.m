function [rot] = get_rotation_from_quaternion(q)
% compute rotation matrix from quaternion
% Input: q = [qW, qX, qY, qZ];

rot = [
    2*(q(1)^2+q(2)^2)-1,     2*(q(2)*q(3)-q(1)*q(4)),    2*(q(2)*q(4)+q(1)*q(3));
    2*(q(2)*q(3)+q(1)*q(4)), 2*(q(1)^2+q(3)^2)-1,        2*(q(3)*q(4)-q(1)*q(2));
    2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)),    2*(q(1)^2+q(4)^2)-1;    
    ];
