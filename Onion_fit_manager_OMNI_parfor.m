close all;
clear all;
load Onion_data_piped(4th1Hz).mat;
format shortg;
startTime = tic;

experiment = 'Onion';
nSubject = 12;
nCondition = 6; % w(0.2 0.6 1) X dv(-0.3 0.3)
% for computing BIC
n = nSubject*(nCondition-1); % like degree of freedom. subject*(condition - 1)
modelList = {
    'nullModel';
    'distanceModel';
    'speedModel';
    'sDistanceModel';
    'expansionModel';
    'expansion2Model';
    'expansion3Model';
    'ratioModel';
    'ratio2Model';
    'ratio3Model';
    'linearModel';
    'lemercierModel'
    };

%                             'nullModel';
%                             'distanceModel';
%                             'speedModel';
%                             'sDistanceModel';
%                             'expansionModel';
%                             'expansion2Model';
%                             'expansion3Model';
%                             'ratioModel';
%                             'ratio2Model';
%                             'ratio3Model';
%                             'linearModel';
%                             'lemercierModel'
                              
% this should correspond to modelList
initParamList = {
    [];
    [1];
    [1];
    [1 1 1];
    [1];
    [1];
    [1 1];
    [1 1 1];
    [1];
    [1 1];
    [1 1 1 1];
    [1 1 1];
    
    };


                        %     [];
                        %     [1];
                        %     [1];
                        %     [1 1 1];
                        %     [1];
                        %     [1];
                        %     [1 1];
                        %     [1 1 1];
                        %     [1];
                        %     [1 1];
                        %     [1 1 1 1];
                        %     [1 1 1];



h = 2;
method = 'ALL'; % MC: Monte Carlo
                % LOO: leave-one-out (leave one trial out of per subjectXcondition)
                % ALL: fit to all data
                % LOSO: leave-one-subject-out
delay = 0;
criterion = 'RMSE';
variable = 'v';
nIter = 1;
if strcmp(method,'LOSO')
    nIter = max([ExperimentalTrials.subject]);
end
simLength = 5.5; % 5.5 seconds of data in each trial is used for training and testing

Hz = 90;
% delete(gcp('nocreate'));
% parpool(4);


%% choose data set
ExperimentalTrials = rmfield(ExperimentalTrials, {'raw_traj';'zero_traj';'rot_traj';'inter_traj';'filtered_traj'});
DATA_SET = ExperimentalTrials(1);

for i = 1:length(ExperimentalTrials)
    if ExperimentalTrials(i).dump ~= 1 && ExperimentalTrials(i).dv ~= 0 &&...
            ExperimentalTrials(i).t_total-ExperimentalTrials(i).manipOnset>simLength
        DATA_SET(end+1) = ExperimentalTrials(i);
    end
end

DATA_SET = DATA_SET(2:end);

%% grouping by participants and conditions if using LOO cross validation
if strcmp(method,'LOO')
    PoolList = struct('subjCon', '', 'trialIDs', [], 'count', 0);
    for subject = 1:nSubject
        for w = [0.2 0.6 1]
            for dv = [-0.3 0.3]
                PoolList(end+1).subjCon = ['subj' num2str(subject) '_' num2str(w) num2str(dv)];
                for i = 1:length(DATA_SET)
                    if DATA_SET(i).subject == subject && DATA_SET(i).w == w &&...
                            DATA_SET(i).dv == dv
                        PoolList(end).trialIDs = [PoolList(end).trialIDs; i];
                    end
                end
            end
            
        end
    end
    PoolList = PoolList(2:end);
    for i = 1:length(PoolList)
        PoolList(i).count = length(PoolList(i).trialIDs);
    end
end

%% fitting and validation
% set up output structure
nTrial = length(DATA_SET);
nModel = length(modelList);
ResultsSummary = struct('model', '', 'iIter', 0, 'initParam', [], 'pIter',...
    [],'fminIter', 0, 'testError', 0, 'iTotalTime', 0);
Results = struct;
Results_ = rmfield(DATA_SET(1),'data');
Results_.tests = NaN;
Results_.z_x = NaN;
Results_.z_v = NaN;
Results_.z_a = NaN;
Results_.RMSE_x = NaN;
Results_.RMSE_v = NaN;
Results_.RMSE_a = NaN;
Results_.model = NaN;
Results_.iIter = NaN;


% generate train/test set split for each iteration
iTrainSet = {};
iTestSet = {};

for i = 1:nIter
    if strcmp(method, 'LOO')
        % evenly draft from each participant and condition
        iTrainSet_ = [];
        iTestSet_ = [];
        for j = 1:length(PoolList)
            iRandTrials = PoolList(j).trialIDs(randperm(PoolList(j).count));
            iTrainSet_ = [iTrainSet_; iRandTrials(2:end)];
            iTestSet_ = [iTestSet_; iRandTrials(1)];        
        end
        iTrainSet(i) = {iTrainSet_};
        iTestSet(i) = {iTestSet_};
    elseif strcmp(method, 'MC')
        iRandTrials = randperm(nTrial);
        iTrainSet(i) = {iRandTrials(1:int16(nTrial*.75))'};
        iTestSet(i) = {iRandTrials(int16(nTrial*.75)+1:end)'};  
    elseif strcmp(method, 'ALL')
        iTrainSet(i) = {(1:nTrial)'};
        iTestSet(i) = {(1:nTrial)'};
    elseif strcmp(method, 'LOSO')
        iTrainSet(i) = {find([DATA_SET.subject]~=i)'};
        iTestSet(i) = {find([DATA_SET.subject]==i)'};
    end
end

indexList = combvec(1:nModel, 1:nIter); 
% frist dimension: the index for model and parameter
% second dimension: iteration
% use the index of one for loop to achieve the effect of the following nested for loop:
%     for iIter = 1:nIter
%         for model = modelList
%             code
%         end
%     end
parfor ii = 1:size(indexList,2)
    
    iStartTime = tic;
    iModel = indexList(1,ii);
    iIter = indexList(2,ii);  
    model = modelList(iModel);
    initParam = initParamList{iModel};
     
    TRAIN_SET = DATA_SET(iTrainSet{iIter});
    TEST_SET = DATA_SET(iTestSet{iIter});

    % train the model
    if ~strcmp(model,'nullModel')
        zFunction = @(p) train(model, p, h, delay, TRAIN_SET, criterion, variable, simLength, Hz);
        [pIter, fminIter] = fminunc(zFunction, initParam);
    else
        pIter = NaN;
        fminIter = NaN;
    end
    
    % test the model and save result to internal variables
    iResults = test(TEST_SET, model, iIter, pIter, h, delay, simLength, Hz);
    ResultsSummary(ii).model = model;
    ResultsSummary(ii).iIter = iIter;
    ResultsSummary(ii).initParam = initParam;
    ResultsSummary(ii).pIter = pIter;
    ResultsSummary(ii).fminIter = fminIter;
    ResultsSummary(ii).testError = mean([iResults.RMSE_v]);
    if ~strcmp(model,'nullModel')
        k = 0;
    else
        k = length(pIter);
    end
    RMSE_v = [iResults.RMSE_v];
    MSE = mean(RMSE_v.^2);
    ResultsSummary(ii).BICWiki = n*log(MSE) + k*log(n);
    Results_ = [Results_, iResults];
    ResultsSummary(ii).iTotalTime = datestr(toc(iStartTime)/(24*60*60), 'HH:MM:SS');
    [char(model), ' on iteration ', num2str(iIter), ' is done.']
end

Results_ = Results_(2:end);
for i=1:length(Results_)
    Results(i).model = Results_(i).model;
    Results(i).iIter = Results_(i).iIter;
end
fn = fieldnames(Results_);
for i = fn(1:end-2)'
    for j = 1:length(Results_)
        Results(j).(i{1}) = Results_(j).(i{1});
    end
end

totalTime = toc(startTime);
totalTime = datestr(totalTime/(24*60*60), 'HH:MM:SS');
time = string(clock);

% write results
save([experiment,'_fit_pertTrials_',method,num2str(nIter),'_',criterion,'_',variable,'_',char(join(time(1:5),'-')),'.mat'],...
     'method','modelList','criterion','variable','delay','nIter','nTrial', 'initParamList',...
     'simLength', 'simLength', 'Results', 'ResultsSummary', 'totalTime');
% pertTrials means perturbation trials dv = -0.3 or 0.3

%% Functions

function index = train(model, p, h, delay, TRAIN_SET, criterion, variable, simLength, Hz)
fisherZs = zeros(length(TRAIN_SET),1);
RMSEs = zeros(length(TRAIN_SET),1);

for i = 1:length(TRAIN_SET)
    manipOnset = TRAIN_SET(i).manipOnset;
    frameStart = int32((manipOnset-0.5)*Hz)+1;
    frameEnd = frameStart + simLength*Hz - 1;
    w = TRAIN_SET(i).w;
    lPos = TRAIN_SET(i).data(frameStart:frameEnd,1);
    lSpd = TRAIN_SET(i).data(frameStart:frameEnd,3);
    fPos = TRAIN_SET(i).data(frameStart:frameEnd,2);
    fSpd = TRAIN_SET(i).data(frameStart:frameEnd,4);
    fAcc = TRAIN_SET(i).data(frameStart:frameEnd,6);
    inputHz = Hz;
    outputHz = Hz;
    pStart = fPos(1);
    vStart = fSpd(1);
    [mPos, mSpd, mAcc] = models(model, p, delay, w, h, lPos, lSpd, pStart, vStart, inputHz, outputHz);
    
    % choose criterion
    if strcmp(criterion,'pearsonR')
        % Use correlation for training
        if strcmp(variable,'x')
            fisherZs(i) = atanh(corr(fPos, mPos, 'type','pearson'));
        end
        if strcmp(variable,'v')
            fisherZs(i) = atanh(corr(fSpd, mSpd, 'type','pearson'));
        end
        if strcmp(variable,'a')
            fisherZs(i) = atanh(corr(fAcc, mAcc, 'type','pearson'));
        end
    elseif strcmp(criterion,'RMSE')
        % Use RMSE for training
        
        if strcmp(variable,'x')
            y = fPos;
            yhat = mPos;
            RMSEs(i) = sqrt(mean((y - yhat).^2));
        end
        if strcmp(variable,'v')
            y = fSpd;
            yhat = mSpd;
            RMSEs(i) = sqrt(mean((y - yhat).^2));

        end
        if strcmp(variable,'a')
            y = fAcc;
            yhat = mAcc;
            RMSEs(i) = sqrt(mean((y - yhat).^2));
        end
    end
end

% choose output
if strcmp(criterion,'pearsonR')
    index = -tanh(mean(fisherZs)); % use negative sign so that the smallest index is the highest correlation
elseif strcmp(criterion,'RMSE')
    index = mean(RMSEs);
end
end


function iResults = test(TEST_SET, model, iIter, p, h, delay, simLength, Hz)
iResults = TEST_SET;
for i = 1:length(iResults)
    manipOnset = iResults(i).manipOnset;
    frameStart = int32((manipOnset-0.5)*Hz)+1;
    frameEnd = frameStart + simLength*Hz - 1;
    w = iResults(i).w;
    lPos = iResults(i).data(frameStart:frameEnd,1);
    lSpd = iResults(i).data(frameStart:frameEnd,3);
    lAcc = iResults(i).data(frameStart:frameEnd,5);
    fPos = iResults(i).data(frameStart:frameEnd,2);
    fSpd = iResults(i).data(frameStart:frameEnd,4);
    fAcc = iResults(i).data(frameStart:frameEnd,6);
    inputHz = Hz;
    outputHz = Hz;
    pStart = fPos(1);
    vStart = fSpd(1);
    
    [mPos, mSpd, mAcc] = models(model, p, delay, w, h, lPos, lSpd, pStart, vStart, inputHz, outputHz);
    
    iResults(i).tests(:,1) = lPos;
    iResults(i).tests(:,2) = fPos;
    iResults(i).tests(:,3) = mPos;
    iResults(i).tests(:,4) = lSpd;
    iResults(i).tests(:,5) = fSpd;
    iResults(i).tests(:,6) = mSpd;
    iResults(i).tests(:,7) = lAcc;
    iResults(i).tests(:,8) = fAcc;
    iResults(i).tests(:,9) = mAcc;
    iResults(i).tests(:,10) = iResults(i).data(frameStart:frameEnd,7);
    iResults(i).z_x = atanh(corr(fPos, mPos, 'type','pearson'));
    iResults(i).z_v = atanh(corr(fSpd, mSpd, 'type','pearson'));
    iResults(i).z_a = atanh(corr(fAcc, mAcc, 'type','pearson'));
    iResults(i).RMSE_x = sqrt(mean((fPos - mPos).^2));
    iResults(i).RMSE_v = sqrt(mean((fSpd - mSpd).^2));
    iResults(i).RMSE_a = sqrt(mean((fAcc - mAcc).^2));
    iResults(i).model = model;
    iResults(i).iIter = iIter;
end

iResults = rmfield(iResults, {'data'}); 

end

function [mPos, mSpd, mAcc] = models(model, p, delay, w, h, lPos, lSpd, pStart, vStart, inputHz, outputHz)
if strcmp(model, 'nullModel')
    [mPos, mSpd, mAcc] = nullModel(pStart, vStart, lSpd, inputHz);
    
elseif strcmp(model,'speedModel')
    c = p(1);
    [mPos, mSpd, mAcc] = speedModel(pStart, vStart, inputHz,outputHz, lSpd, c, delay);
    
elseif strcmp(model, 'distanceModel')
    c = p(1);
    [mPos, mSpd, mAcc] = distanceModel(pStart, vStart, inputHz,outputHz, lPos, c, delay);

    
elseif strcmp(model, 'sDistanceModel')
    c = p(1);
    alpha = p(2);
    beta = p(3);
    [mPos, mSpd, mAcc] = sDistanceModel(pStart, vStart, inputHz,outputHz, lPos, c, alpha, beta, delay);
    
elseif strcmp(model, 'ratioModel')
    c = p(1);
    M = p(2);
    L = p(3);
    [mPos, mSpd, mAcc] = ratioModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, c, M, L, delay);
    
elseif strcmp(model, 'linearModel')
    c1 = p(1);
    c2 = p(2);
    alpha = p(3);
    beta = p(4);
    [mPos, mSpd, mAcc] = linearModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, c1, c2, alpha, beta, delay);
    
elseif strcmp(model, 'lemercierModel')
    RT = p(1);
    C = p(2);
    gamma = p(3);
    [mPos, mSpd, mAcc] = lemercierModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, RT, C, gamma);
    
elseif strcmp(model, 'bruneauModel')
    ttr = p(1);
    df = p(2);
    [mPos, mSpd, mAcc] = bruneauModel(pStart, vStart,inputHz,outputHz, lPos, ttr, df);
    
elseif strcmp(model, 'ratio2Model')
    C = p(1);
    [mPos, mSpd, mAcc] = ratio2Model(pStart, vStart,inputHz,outputHz, lPos, lSpd, C, delay);
    
elseif strcmp(model, 'ratio3Model')
    C = p(1);
    L = p(2);
    [mPos, mSpd, mAcc] = ratio3Model(pStart, vStart,inputHz,outputHz, lPos, lSpd, C, L, delay);
         
elseif strcmp(model, 'expansionModel')
    b = p(1);
    [mPos, mSpd, mAcc] = expansionModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay);
    
elseif strcmp(model, 'expansion2Model')
    b = p(1);
    [mPos, mSpd, mAcc] = expansion2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay);
    
elseif strcmp(model, 'expansion3Model')
    b = p(1);
    L = p(2);
    [mPos, mSpd, mAcc] = expansion3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, L, delay);
    
elseif strcmp(model, 'expansionHModel')
    b = p(1);
    [mPos, mSpd, mAcc] = expansionHModel( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, delay);
    
elseif strcmp(model, 'expansionH2Model')
    b = p(1);
    [mPos, mSpd, mAcc] = expansionH2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, delay);
    
elseif strcmp(model, 'expansionH3Model')
    b = p(1);
    L = p(2);
    [mPos, mSpd, mAcc] = expansionH3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, L, delay);
    
elseif strcmp(model, 'expansionWHModel')
    b = p(1);
    [mPos, mSpd, mAcc] = expansionWHModel( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, h, delay);
    
elseif strcmp(model, 'expansionWH2Model')
    b = p(1);
    [mPos, mSpd, mAcc] = expansionWH2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, h, delay);
    
elseif strcmp(model, 'expansionWH3Model')
    b = p(1);
    L = p(2);
    [mPos, mSpd, mAcc] = expansionWH3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, h, L, delay);
    
end
end

function [fPos, fSpd, fAcc] = nullModel(pStart, vStart, lSpd, inputHz)
% This function simulate a follower with no acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                                   x..f = 0
% 
% x..f: the current acceleration of follower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(lSpd);

fAcc = zeros(n,1);
fSpd = repmat(vStart, n, 1);
fPos = pStart:fSpd/inputHz:(pStart + (n-1)*fSpd/inputHz);
fPos = fPos';


end


function [fPos, fSpd, fAcc] = speedModel(pStart, vStart,inputHz,outputHz, lSpd, c, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the speed change of the
% leader.
% c: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                           x..f = c * [x.l - x.f]
% 
% x..f: current acceleration of follower
% c: free parameter
% x.l: current speed of leader
% x.f: current speed of follower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lSpd = ExperimentalTrials(1).data(:,3);
% c = 1.6;
% pStart = 0;
% vStart = 0;
% inputHz = 60;
% outputHz = 60;
Hz = 6000;
time = length(lSpd)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = distanceModel(pStart, vStart,inputHz,outputHz, lPos, c, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% c: free parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                     x..f = c * [delta_x - delta_x0]
% 
% x..f: the current acceleration of follower
% c: free parameter
% delta_x: the current distance between leader and follower
% delta_x0 = d0: initial distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!    
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lPosExtrap(i-int32(delay*Hz)) - fPos(i) - lPosExtrap(1));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = sDistanceModel (pStart, vStart,inputHz,outputHz, lPos, c, alpha, beta, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the speed change of the
% leader.
% c, alpha, beta: free parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                  x..f = c * [delta_x - alpha - beta*x.f]
% 
% x..f: the current acceleration of follower
% c, alpha, beta: free parameter
% delta_x: the current distance between leader and follower
% x.f: the current speed of the follower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if beta < 0
    beta = 0;
end
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * (lPosExtrap(i-int32(delay*Hz)) - fPos(i) - alpha - beta * fSpd(i));
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = ratioModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, c, M, L, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% c, M, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f(t) = c * x.f^M * delta_x. / delta_x^L
% 
% x..f(t): the current acceleration of follower
% x.f: the speed of the follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% c, M, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c * fSpd(i)^M * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i)) / (lPosExtrap(i-int32(delay*Hz)) - fPos(i))^L; 
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = linearModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, c1, c2, alpha, beta, delay)
% This function can simulate the trajectory of follower given that of the
% leader. 
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% c1, c2, alpha, beta: free parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%            x..f = c1*[delta_x.] + c2*[delta_x - alpha - beta*x.f]
% 
% x..f: current acceleration of follower
% delta_x.: relative speed
% delta_x: relative distance
% x.f: follower speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if beta < 0
    beta = 0;
end
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = c1 * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i)) + c2 * (lPosExtrap(i-int32(delay*Hz)) - fPos(i) - alpha - beta * fSpd(i));
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = lemercierModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, RT, C, gamma)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% RT: the delay before the follower react to the distance change between
% the follower and the leader.
% C, gamma, RT: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x..f(t) = C * delta_x.(t - RT) * (1/delta_x(t))^gamma
% 
% x..f(t): the current acceleration of follower
% x.f: the speed of the follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% C, gamma: free parameter
% RT: reaction time, range [0 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RT < 0
    RT = 0;
end
if RT > 1
    RT = 1;
end

Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(RT*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = C * (lSpdExtrap(i-int32(RT*Hz)) - fSpd(i-int32(RT*Hz))) * (1/(lPosExtrap(i) - fPos(i)))^gamma;
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end



function [fPos, fSpd, fAcc] = bruneauModel(pStart, vStart,inputHz,outputHz, lPos, ttr, df)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x.f(t) = (x_l(t+delta_t) - x_f(t) - df) / (delta_t + ttr)
% 
% x.f: follower's speed (m/s)
% delta_t: time step (s)
% ttr: time to react (s)
% df: sum of body size and personal distance (m)
% x_f: follower's position (m)
% x_l: leader's position (m)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ttr < 0
    ttr = 0;
end
if ttr > 1
    ttr = 1;
end
timeStep = 0.1; % second

Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2:(nFrame-timeStep*Hz)
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i-1) + fSpd(i-1)/Hz;
    fSpd(i) = (lPosExtrap(i+ timeStep*Hz) - fPos(i) - df)/(timeStep + ttr);    
    fAcc(i) = (fSpd(i) - fSpd(i-1))*Hz;
end

fPos((nFrame-timeStep*Hz+1):nFrame) = repmat(fPos(nFrame-timeStep*Hz), timeStep*Hz, 1);
fSpd((nFrame-timeStep*Hz+1):nFrame) = repmat(fSpd(nFrame-timeStep*Hz), timeStep*Hz, 1);
fAcc((nFrame-timeStep*Hz+1):nFrame) = repmat(fAcc(nFrame-timeStep*Hz), timeStep*Hz, 1);
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;

% fPos = fPos(1:100:6000*6);
% fSpd = fSpd(1:100:6000*6);
% fAcc = fAcc(1:100:6000*6);

fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';



end

function [fPos, fSpd, fAcc] = ratio2Model(pStart, vStart,inputHz,outputHz, lPos, lSpd, C, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x..f(t) = C * delta_x. / delta_x
%
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% x..f(t): the current acceleration of follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% C, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = C * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i-int32(delay*Hz))) / (lPosExtrap(i) - fPos(i));
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = ratio3Model(pStart, vStart,inputHz,outputHz, lPos, lSpd, C, L, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%          x..f(t) = C * delta_x. / delta_x^L
%
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% x..f(t): the current acceleration of follower
% delta_x: the current distance between leader and follower
% delta_x.: the rate of change of the current distance between leader and follower
% C, L: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');
% because this study is about locomotion on a line so only 1 dimension is considered.
fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = C * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i-int32(delay*Hz))) / (lPosExtrap(i) - fPos(i))^L;
end
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = expansionModel( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha.
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target width
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * w * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + w^2/4);    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end

function [fPos, fSpd, fAcc] = expansion2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha./alpha
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target width
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * w * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/...
            ((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + w^2/4)/...
            (2*atan(w/2/(lPosExtrap(i-int32(delay*Hz)) - fPos(i))));     
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end



function [fPos, fSpd, fAcc] = expansion3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, L, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha./alpha^L
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target width
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * w * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/...
            ((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + w^2/4)/...
            (2*atan(w/2/(lPosExtrap(i-int32(delay*Hz)) - fPos(i))))^L;    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = expansionHModel( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha.
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target height
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * h * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + h^2/4);    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = expansionH2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha./alpha
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target height
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * h * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/...
            ((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + h^2/4)/...
            (2*atan(h/2/(lPosExtrap(i-int32(delay*Hz)) - fPos(i))));     
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = expansionH3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, h, L, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * alpha./alpha^L
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target height
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * h * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/...
            ((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + h^2/4)/...
            (2*atan(h/2/(lPosExtrap(i-int32(delay*Hz)) - fPos(i))))^L;    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end


function [fPos, fSpd, fAcc] = expansionWHModel( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, h, delay)
% This function can simulate the trajectory of follower given that of the
% leader.
% f: follower
% l: leader
% pos: position   spd: speed    acc: acceleration
% b: free parameter
% w, h: width and height of the target
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The model it use:  
% 
%                x..f = -b * [alpha(w)*alpha(h)].
% 
% x..f(t): the current acceleration of follower
% alpha.: the changing rate of visual angle based on target width * 
% b: free parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hz = 6000;
time = length(lPos)/inputHz;
nFrame = int32(time*Hz);
x = 1/inputHz:1/inputHz:time;
xq = 1/Hz:1/Hz:time;
lPosExtrap = interp1(x,lPos,xq,'linear','extrap');
lSpdExtrap = interp1(x,lSpd,xq,'linear','extrap');

fPos = repmat(pStart, nFrame, 1);
fSpd = repmat(vStart, nFrame, 1);
fAcc = zeros(nFrame,1);

for i = 2+int32(delay*Hz):nFrame
    % the units in subscript should be frame
    % Be sure to update each variable by updated variable!!!!!!!!!!!
    % Is Acc is updated by the current fSpd and fPos, it should be the last
    % line!
    fPos(i) = fPos(i - 1) + fSpd(i - 1)/Hz;
    fSpd(i) = fSpd(i - 1) + fAcc(i - 1)/Hz;
    fAcc(i) = b * h * (lSpdExtrap(i-int32(delay*Hz)) - fSpd(i))/((lPosExtrap(i-int32(delay*Hz)) - fPos(i))^2 + h^2/4);    
end 
x = 1/Hz:1/Hz:time;
xq = 1/outputHz:1/outputHz:time;
fPos = interp1(x,fPos,xq,'linear','extrap')';
fSpd = interp1(x,fSpd,xq,'linear','extrap')';
fAcc = interp1(x,fAcc,xq,'linear','extrap')';
end