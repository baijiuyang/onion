
nDuration = 12; % 12 second, the duration of the trajectory
startupDuration = 0; % how fast the pole start to move
meanManipOnset = 2.5; % mean manipulation oneset time
onsetWindow = 1; % the window of manipulation onset

frameRate = 90;

% part of the input to trajectoryGenerator
heading1 = 0; % pre-manipulation heading
heading2 = 0; % post-manipulation heading
x0 = 0; % the pole is moving on y dimension so x0 is constant
d0 = 4;
v0 = 0.8;
dv = 0.3;
a = 1;

[x, y, spd, hdn, manipOnset] = trajectoryGenerator(x0,d0,nDuration,v0,dv,a,...
    heading1,heading2,startupDuration,meanManipOnset,onsetWindow,frameRate);

speed = diff(y)*frameRate;
figure;
plot(1:length(speed),speed);
csvwrite('test.csv', y);

%% plot test csv
figure;
file = load('test.csv');
traj = file(:,1);
speed = diff(traj)*90;
plot(1:length(speed),speed);

%% plot loaded trial
figure;
for i = 1
    file = load(strcat('trial',num2str(i,'%03d'),'.csv'));
    traj = file(:,2);
    speed = diff(traj)*90;
    hold on;
    plot(1:length(speed),speed);
end
