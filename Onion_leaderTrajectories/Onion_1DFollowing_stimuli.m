% =======================================================================
% GENERATE EXPERIMENT J Bai Experiment; Speed Manipulation Pole following
%       Created from Kevin Rio's original Density Experiment code
%       Altered by TW JB
%
% Created on: January 02 2017
% Revised on: April 29 2017 simplify the loop structure
% Revised on: May 03 2017 add manipOnset for output
% Revised on: 2/19/2018 by JB. Write all input file in one folder
% instead of separate folders. Generate all subjects' stimuli by one click.
% =======================================================================

close all;
% parameters and initialization
totalSubject = 15;
EXPERIMENT = 'Onion';
for subjectID = 0:totalSubject

    nRepetition = 10;
    nDuration = 12; % 12 second, the duration of the trajectory
    startupDuration = 0; % how fast the pole start to move
    meanManipOnset = 3.5; % mean manipulation oneset time
    onsetWindow = 1; % the window of manipulation onset (second)

    frameRate = 90;

    % part of the input to trajectoryGenerator
    heading1 = 0; % pre-manipulation heading
    heading2 = 0; % post-manipulation heading
    x0 = 0; % the pole is moving on y dimension so x0 is constant
    % variables and levels
    d0 = [2];
    v0 = [1.2];
    dv = [-0.3, 0, 0.3];
    a = [1];
    size = [0.2 0.6 1];
    % number of levels in each variable
    nSize = length(size);
    nDv = length(dv);


    iTrial = 1; % counter of simualtions
    nTrial = nSize*nDv*nRepetition; % total number of trials

    % data of each frame
    Nfo = cell(nTrial,1); % trial information
    Pos = cell(nTrial,1); % pole position
    Spd = cell(nTrial,1); % pole speed
    Hdn = cell(nTrial,1); % pole heading

    % The nested simulation loops
    for i = 1:nSize
        for j = 1:nDv
            for k = 1:nRepetition

                Nfo{iTrial}.d0 = d0;
                Nfo{iTrial}.v0 = v0;
                Nfo{iTrial}.dv = dv(j);
                Nfo{iTrial}.a = a;


                [x,y,spd,hdn, manipOnset] = inputTrajectoryGenerator(x0,d0,nDuration,v0,dv(j),a,...
                    heading1,heading2,startupDuration,meanManipOnset,onsetWindow,frameRate);


                Nfo{iTrial}.manipOnset = manipOnset; % the start time of speed change of target pole
                Nfo{iTrial}.size = size(i);
                Pos{iTrial,1}(:,1) = x;
                Pos{iTrial,1}(:,2) = y;
                Spd{iTrial,1} = spd;
                Hdn{iTrial,1} = hdn;

                % update the counter of simulations
                if iTrial < nTrial
                    iTrial = iTrial + 1;
                end
            end
        end
    end


    % shuffle the trials by the 'SuffleTrials' function
    [shuffledNfo,shuffledPos,shuffledSpd,shuffledHdn] = shuffleTrials(Nfo, Pos, Spd, Hdn);



    % wrtie the CSV file
    writeTrajectory(EXPERIMENT, subjectID, shuffledNfo, shuffledPos, shuffledSpd, frameRate);

    % unshuffled trials use the following code:
    % writeTrajectory(Nfo, Pos, Spd, frameRate);

end
