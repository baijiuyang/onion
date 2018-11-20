function [x, y, spd, hdn, manipOnset] = inputTrajectoryGenerator(x0,y0,nDuration,v0,dv,a,...
    heading1,heading2,startupDuration,meanManipOnset,onsetWindow,frameRate)
% VIRTUAL PEDESTRIAN TRAJECTORY Generates virtual pedestrian trajectory
%   Used for Kevin's Virtual Crowd experiment. Given a large set of initial
%   conditions and parameters, generates the time series of x, y, speed,
%   and heading. Used to pre-program the trajectories of virtual
%   pedestrians in a large crowd. 
%
%   Trajectories are divided into 4 phases: initial acceleration
%   (transition from standstill to v0), pre-manipulation (steady-state 
%   at v0), manipulation (transition from v0 to v0 + dv), and
%   post-manipulation (steady-state at v0 + dv). Transitions are smoothed
%   using a normal cumulative density function that has been calibrated to
%   human data. 
%
%   INPUTS
%       x0:              initial x-position (m)
%       y0:              initial y-position (m)
%       nDuration:           number of time steps to generate (1/frameRate s)
%       v0:          speed pre-manipulation (m/s)
%       dv:          delta_v speed change (m/s)
%       heading1:        heading pre-manipulation (deg)
%       heading2:        heading post-manipulation (deg)
%       startupDuration: length of initial acceleration ("startup") period (1/frameRate s)
%       manipTime:       time at which manipulation begins (1/frameRate s)
%       manipDuration:   length of manipulation period (1/frameRate s)
%       frameRate:       how many frames in a second
%
%   OUTPUTS
%       x:               x-coordinate of position (m) vector
%       y:               y-coordinate of position (m) vector
%       spd:             speed time series (m/s) vector
%       hdn:             heading time series (deg) vector

% Compute normal cumulative density function, to smoothly transition
% between pre- and post-manipulation speed/heading. 

% nTime was changed to nDuration_Jiuyang Bai
% manipTime is randomized around a mean_Jiuyang Bai


manipDuration = abs(dv) / a;

% convert unit from (s), (m/s) to (frame), (m/frame) 
v0 = v0 * (1/frameRate);
dv = dv * (1/frameRate);
nDuration = nDuration * frameRate;
if startupDuration == 0
    startupDuration = 1;
else
    startupDuration = startupDuration * frameRate;
end
meanManipOnset = meanManipOnset * frameRate;
manipDuration = manipDuration * frameRate;
% random start time of manipulation (+- onsetWindow/2 s) around a mean
manipOnset = meanManipOnset + int16((rand - 0.5)*onsetWindow*frameRate);
heading1 = heading1 * (pi/180);
heading2 = heading2 * (pi/180);


% gradual acceleration
startup = normcdf(-startupDuration/2:startupDuration/2,0,startupDuration/6);
manip = normcdf(-manipDuration/2:manipDuration/2,0,manipDuration/6);


% Initial conditions.
x = NaN(nDuration,1);
y = NaN(nDuration,1);
x(1) = x0;
y(1) = y0;
hdn = NaN(nDuration,1);
spd = NaN(nDuration,1);
hdn(1) = heading1;


% Main loop.
for iTime = 2:nDuration
    
    % Initial acceleration.
    if iTime <= startupDuration
        hdn(iTime) = heading1;
        x(iTime) = x(iTime-1);
        y(iTime) = y(iTime-1) + v0*startup(iTime);
        
    % Pre-manipulation.
    elseif iTime > startupDuration && iTime <= manipOnset
        hdn(iTime) = heading1;
        x(iTime) = x(iTime-1);
        y(iTime) = y(iTime-1) + v0;
        
    % Manipulation.
    elseif iTime > manipOnset && iTime <= manipOnset + manipDuration
        hdn(iTime) = (heading2-heading1)*manip(iTime-manipOnset) + heading1;
        x(iTime) = x(iTime-1) + ... 
            (dv*manip(iTime-manipOnset) + v0) ...
            *sin(hdn(iTime));
        y(iTime) = y(iTime-1) + ... 
            (dv*manip(iTime-manipOnset) + v0) ...
            *cos(hdn(iTime));            
        
    % Post-manipulation.
    elseif iTime > manipOnset + manipDuration
        hdn(iTime) = heading2;
        x(iTime) = x(iTime-1) + (v0 + dv)*sin(heading2);
        y(iTime) = y(iTime-1) + (v0 + dv)*cos(heading2);     
    end
 
end

hdn = hdn*(180/pi);
spd = sqrt(diff(x).^2+diff(y).^2)*(frameRate);
spd(nDuration) = spd(nDuration-1);


end