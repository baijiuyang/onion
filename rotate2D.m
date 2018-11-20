function new_traj = rotate2D(theta, traj, origin)
	% rotation function that rotates a trajectory around a given origin point
    % theta is in degrees
    % Jiuyang Bai Sep. 1st 2017
    
    new_traj = traj - origin;
    
    R = [cosd(theta), sind(theta);...
         -sind(theta), cosd(theta)];
    
    new_traj = new_traj * R;
    
    new_traj = new_traj + origin;
    
