function [kspLocOut, phaseOut] = traj_correction_GIRF(kspLocIn, dt, rot, acqDelay, gradDelay)

% gamma = 4257.59;
gamma = 4257.7;
% gamma = 2.67522212e+8;
R.T = 0.55;
R.R = rot;

grad = cat(2,kspLocIn(:,1,:), diff(kspLocIn,1,2))/dt/(gamma);
grad = permute(grad,[2 3 1]);
grad(:,:,3) = 0;

% [kspLocOut, ~] = apply_GIRF(grad, dt, R, delay_offset); phaseOut = [];
[kspLocOut, ~, phaseOut] = apply_GIRF_20210827(grad, dt, R, acqDelay, gradDelay);
kspLocOut = permute(kspLocOut, [3, 1, 2]);

