clear all;
leader_traj_source = 'perfect'; % 'perfect'/'vizard_input'/'vizard_output'

cutoff = 1; % cutoff frequency of the butterworth filter
order = 1; % the order of the cutoff frequency
pad = 180; % length of data extension before filtering

%% loading
subFolder = 'Onion_rawData\';

% load test trials
loading1 = dir(strcat(subFolder,'Onion_subj*.csv')); % load all txt files in directory of code
ExperimentalTrials = struct;


for iTrial = 1:length(loading1) % iterate through all trials
    % read subject number
    s_ = strfind(loading1(iTrial).name, 'subj');
    ExperimentalTrials(iTrial).subject = str2double(loading1(iTrial).name(s_+4 : s_+5));
    
    % read trial number
    s_ = strfind(loading1(iTrial).name, 'trial');
    ExperimentalTrials(iTrial).trial = str2double(loading1(iTrial).name(s_+5 : s_+7));
    
    % read condition
    ExperimentalTrials(iTrial).w = str2double(loading1(iTrial).name(23:25));
    ExperimentalTrials(iTrial).dv = str2double(loading1(iTrial).name(end-7:end-4));
    if isnan(ExperimentalTrials(iTrial).dv)
        s_ = strfind(loading1(iTrial).name, '_DUMP_');
        ExperimentalTrials(iTrial).dv = str2double(loading1(iTrial).name(s_-4:s_-1));
    end

    fileName_ = strcat(subFolder, 'Onion_input\Onion_subject', num2str(ExperimentalTrials(iTrial).subject,'%02d'),...
        '_conditions.csv');
    conditionFile = csvread(fileName_,1,0); % read data start from row 2 column 1
    % make sure each trial is reading the right manipulation onset time by
    % comparing the w and dv of that trial
    if ExperimentalTrials(iTrial).w == conditionFile(ExperimentalTrials(iTrial).trial,6) &&...
       ExperimentalTrials(iTrial).dv == conditionFile(ExperimentalTrials(iTrial).trial,4)
        ExperimentalTrials(iTrial).manipOnset = conditionFile(ExperimentalTrials(iTrial).trial,5);
    end

    fileName_ = strcat(subFolder, loading1(iTrial).name); % find name of i_trial file
    raw_traj = load(fileName_);
    ExperimentalTrials(iTrial).t_total = raw_traj(end,end);
    ExperimentalTrials(iTrial).raw_traj = [raw_traj(:,1:4) raw_traj(:,end)];
    
end

% % load freewalk trials
% loading2 = dir(strcat(subFolder,'position_freewalk*')); % load all txt files in directory of code
% freewalk = struct;
% 
% for iTrial = 1:length(loading2) % iterate through all trials
%     
%     s_ = strfind(loading2(iTrial).name, 'subj');
%     freewalk(iTrial).subject = str2double(loading2(iTrial).name(s_+5 : s_+6));
% 
%     freewalk(iTrial).session = loading2(iTrial).name(20);
%     
%     if strcmp(loading2(iTrial).name(end-5), '_')
%         freewalk(iTrial).trial = str2double(loading2(iTrial).name(end-4));
%     else
%         freewalk(iTrial).trial = str2double(loading2(iTrial).name(end-5:end-4));
%     end
%     
%     fileName_ = strcat(subFolder, loading2(iTrial).name); % find name of i_trial file
%     freewalk(iTrial).raw_traj = load(fileName_); % load in specific file
%     
% end


% load error notes
loading3 = dir(strcat(subFolder,'*DUMP*')); % load all txt files in directory of code
Dump = struct;

for iTrial = 1:length(loading3)
    s_ = strfind(loading3(iTrial).name, 'subj');
    Dump(iTrial).subject = str2double(loading3(iTrial).name(s_+4:s_+5));
    s_ = strfind(loading3(iTrial).name, 'trial');
    Dump(iTrial).trial = str2double(loading3(iTrial).name(s_+5:s_+7));
    s_ = strfind(loading3(iTrial).name, '_DUMP_');
    Dump(iTrial).error = loading3(iTrial).name(s_+6:end-4);   
end


%% processing following data
theta = atand(9/11);
Hz = 90;
t = 2; % used the last 2 second to calculate ending speed
sample = Hz*t; % sample of data points used to compute ending speed
v0 = 1.2;
d0 = 2;
a = 1;

for iTrial = 1:length(ExperimentalTrials)
    
    subject = ExperimentalTrials(iTrial).subject;
    trial = ExperimentalTrials(iTrial).trial;
    w = ExperimentalTrials(iTrial).w;
    dv = ExperimentalTrials(iTrial).dv;
    manipOnset = ExperimentalTrials(iTrial).manipOnset;
    
    % make initial position zero
    pos_0 = ExperimentalTrials(iTrial).raw_traj(1,3:4);
    ExperimentalTrials(iTrial).zero_traj(:, 1:2) = ExperimentalTrials(iTrial).raw_traj(:, 1:2) - pos_0;
    ExperimentalTrials(iTrial).zero_traj(:, 3:4) = ExperimentalTrials(iTrial).raw_traj(:, 3:4) - pos_0;
    ExperimentalTrials(iTrial).zero_traj(:, 5) = ExperimentalTrials(iTrial).raw_traj(:, 5); % copy time stamp
    

    % rotate the trajectory to 1 dimension (y)
    ExperimentalTrials(iTrial).rot_traj(:, 1:2) = rotate2D(theta, ExperimentalTrials(iTrial).zero_traj(:, 1:2), [0, 0]);
    ExperimentalTrials(iTrial).rot_traj(:, 3:4) = rotate2D(theta, ExperimentalTrials(iTrial).zero_traj(:, 3:4), [0, 0]);
    % unify two directions of walking
    if ExperimentalTrials(iTrial).rot_traj(end,4) < 0
       ExperimentalTrials(iTrial).rot_traj(:, 1:4) = ExperimentalTrials(iTrial).rot_traj(:, 1:4)*(-1); 
    end
    ExperimentalTrials(iTrial).rot_traj(:, 5) = ExperimentalTrials(iTrial).raw_traj(:, 5); % copy time stamp
    
    % interpolate the trajectory to 90 Hz data with equal time interval 
    x = ExperimentalTrials(iTrial).raw_traj(:, 5);
    v = ExperimentalTrials(iTrial).rot_traj(:, 1:4);
    xp = 1/Hz:1/Hz:ExperimentalTrials(iTrial).raw_traj(end,5);
    vq = interp1(x,v,xp,'linear','extrap');
    ExperimentalTrials(iTrial).inter_traj(:, 1:4) = vq;
    ExperimentalTrials(iTrial).inter_traj(:, 5) = xp;
 
    if strcmp(leader_traj_source, 'vizard_input')
        % load and use leader's trajectory from input files
        fileName_ = strcat(subFolder, 'Carrot4_input\Carrot4_subject', num2str(subject,'%02d'),...
            '_trial',num2str(trial,'%03d'),'.csv');
        inputFile = load(fileName_);
        input_traj = inputFile(:,2);
        trial_length = 12;
        Hz = 90;

        % replace leader's trajectory from vizard by that from the input file
        len = size(ExperimentalTrials(iTrial).inter_traj,1);
        ExperimentalTrials(iTrial).inter_traj(:,2) = input_traj(1:len);
   
    end
    
    
    if strcmp(leader_traj_source, 'perfect')
        
        x0 = 0;
        y0 = 0;
        nDuration = 13;
        heading1 = 0;
        heading2 = 0;
        startupDuration = 0; % how fast the pole start to move
        [x, y, spd, hdn] = Carrot4_trajectoryGenerator(x0,y0,nDuration,...
            d0,v0,dv,a,heading1,heading2,startupDuration,manipOnset, Hz);
        len_ = size(ExperimentalTrials(iTrial).inter_traj,1);
        ExperimentalTrials(iTrial).inter_traj(:,2) = y(1:len_);
    end
    
    % Butterworth filter
    [vOutput, outFiltered, outExtended] = Onion_filter_butter (Hz,ExperimentalTrials(iTrial).inter_traj(:, 4),cutoff,order,pad);
    ExperimentalTrials(iTrial).filtered_traj(:, 2) = vOutput;
    ExperimentalTrials(iTrial).filtered_traj(:, 1) = ExperimentalTrials(iTrial).inter_traj(:, 2); % copy leader's traj
    ExperimentalTrials(iTrial).filtered_traj(:, 3) = ExperimentalTrials(iTrial).inter_traj(:, 5); % copy time stamp after interpolation
    
    
    % calculate speed and acceleration
    % position
    ExperimentalTrials(iTrial).data(:, 1) = ExperimentalTrials(iTrial).filtered_traj(:, 1);% leader
    ExperimentalTrials(iTrial).data(:, 2) = ExperimentalTrials(iTrial).filtered_traj(:, 2);% follower
    
    % speed
    d = diff(ExperimentalTrials(iTrial).data(:, 1));
    ExperimentalTrials(iTrial).data(:, 3) = [d; d(end)] * Hz;% leader
    d = diff(ExperimentalTrials(iTrial).data(:, 2));
    ExperimentalTrials(iTrial).data(:, 4) = [d; d(end)] * Hz;% follower
    
    % acceleration
    d = diff(ExperimentalTrials(iTrial).data(:, 3));
    ExperimentalTrials(iTrial).data(:, 5) = [d; d(end)] * Hz;% leader
    d = diff(ExperimentalTrials(iTrial).data(:, 4));
    ExperimentalTrials(iTrial).data(:, 6) = [d; d(end)] * Hz;% follower
    
    % time stamp
    ExperimentalTrials(iTrial).data(:, 7) = ExperimentalTrials(iTrial).filtered_traj(:, 3);
    
    
    % mark thrown trials
    % mark all trials as good first
    ExperimentalTrials(iTrial).dump = 0;
    
    % mark dump trials by error notes
    for j_trial = 1:length(Dump)
        if ExperimentalTrials(iTrial).subject == Dump(j_trial).subject && ExperimentalTrials(iTrial).trial == Dump(j_trial).trial
            ExperimentalTrials(iTrial).dump = 1;  
        end
    end
    
    % mark dump trials by suspicious speed in filtered data
    if ExperimentalTrials(iTrial).dump == 0 && sum(ExperimentalTrials(iTrial).data(4*Hz:end-Hz,4)<0.1)>0
        ExperimentalTrials(iTrial).dump = 2;
        Dump(end+1).subject = subject;
        Dump(end).trial = trial;
        Dump(end).error = 'suspicious speed';     
    end
    
    % mark dump trials by big speed fluctuation in raw data
    if ExperimentalTrials(iTrial).dump == 0
        speed = diff(ExperimentalTrials(iTrial).inter_traj(:,4));
        speed = [speed; speed(end)];
        speed = speed * 90;
    end
    if ExperimentalTrials(iTrial).dump == 0 && (sum(speed<-5)>0 ||...
            sum(speed>20)>0)
        ExperimentalTrials(iTrial).dump = 2;
        Dump(end+1).subject = subject;
        Dump(end).trial = trial;
        Dump(end).error = 'abnormal fluctuation';     
    end
    
    % ending speed
    ExperimentalTrials(iTrial).finalSpd = mean(ExperimentalTrials(iTrial).data(end-sample+1:end, 4));
    ExperimentalTrials(iTrial).finalDv = mean(ExperimentalTrials(iTrial).data(end-sample+1:end, 3)...
        - ExperimentalTrials(iTrial).data(end-sample+1:end,4));
    ExperimentalTrials(iTrial).finalDist = mean(ExperimentalTrials(iTrial).data(end-sample+1:end,1)...
        - ExperimentalTrials(iTrial).data(end-sample+1:end,2));

end

%% Processing freewalk

% for iTrial = 1:length(freewalk)
%     
%     % make initial position zero
%     pos_0 = freewalk(iTrial).raw_traj(1,3:4);
%     freewalk(iTrial).zero_traj(:, 1:2) = freewalk(iTrial).raw_traj(:, 3:4) - pos_0;
%     freewalk(iTrial).zero_traj(:, 3) = freewalk(iTrial).raw_traj(:, 5); % copy time stamp
%     
%     % rotate the trajectory to 1 dimension (y)
%     freewalk(iTrial).rot_traj(:, 1:2) = rotate2D(theta, freewalk(iTrial).zero_traj(:, 1:2), [0, 0]);
%     if freewalk(iTrial).rot_traj(end,2) < 0
%        freewalk(iTrial).rot_traj(:, 1:2) = freewalk(iTrial).rot_traj(:, 1:2)*(-1); 
%     end
%     freewalk(iTrial).rot_traj(:, 3) = freewalk(iTrial).raw_traj(:, 5); % copy time stamp
%     
%     % interpolate the trajectory to 60 Hz data with equal time interval 
%     x = freewalk(iTrial).raw_traj(:, 5);
%     v = freewalk(iTrial).rot_traj(:, 1:2);
%     xp = 1/Hz:1/Hz:freewalk(iTrial).raw_traj(end,5);
%     vq = interp1(x,v,xp,'linear','extrap');
%     freewalk(iTrial).inter_traj(:, 1:2) = vq;
%     freewalk(iTrial).inter_traj(:, 3) = xp;
%     
%     % Butterworth filter
%     freewalk(iTrial).filtered_traj(:, 1) = filter_butter(Hz, freewalk(iTrial).inter_traj(:, 2));
%     freewalk(iTrial).filtered_traj(:, 2) = freewalk(iTrial).inter_traj(:, 3); % copy time stamp after interpolation
%     
%     
%     % calculate speed and acceleration
%     % position
%     freewalk(iTrial).data(:, 1) = freewalk(iTrial).filtered_traj(:, 1);
%     
%     % speed
%     d = diff(freewalk(iTrial).data(:, 1));
%     freewalk(iTrial).data(:, 2) = [d; d(end)] * Hz;
%     
%     % acceleration
%     d = diff(freewalk(iTrial).data(:, 2));
%     freewalk(iTrial).data(:, 3) = [d; d(end)] * Hz;
%     
%     % time stamp
%     freewalk(iTrial).data(:, 4) = freewalk(iTrial).filtered_traj(:, 2);
%     
% end



%% save
time = string(clock);
save(['Onion_data_piped(', num2str(order), 'th', num2str(cutoff), 'Hz)',char(join(time(1:5),'-')) ,'.mat'], 'ExperimentalTrials', 'Dump');



