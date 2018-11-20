
clear all;
load Carrot4_fit_pertTrials_ALL1_RMSE_v_2018-4-24-17-20(4th1cut).mat;

%% set up the structures


Hz = 90;
time = 5.5; %seconds
nFrame = time*Hz;


ResultsAveraged = struct;
for i = 1:length(modelList)
    for j = 0:12
        for k = [0.2 0.6 1]
            for l = [-0.3 0.3]
                ResultsAveraged(end+1).model = modelList(i);
                ResultsAveraged(end).subject = j;
                ResultsAveraged(end).w = k;
                ResultsAveraged(end).dv = l;
                ResultsAveraged(end).tests = zeros(nFrame,10);
                ResultsAveraged(end).fPos = [];
                ResultsAveraged(end).fSpd = [];
                ResultsAveraged(end).fAcc = [];
                ResultsAveraged(end).mPos = [];
                ResultsAveraged(end).mSpd = [];
                ResultsAveraged(end).mAcc = [];
                ResultsAveraged(end).n = 0;
            end
        end
    end
end

ResultsAveraged = ResultsAveraged(2:end);


%% Calculate average data and tests and compute confidence interval

% loop through all trials to get the sum
for i = 1:length(Results)
    model = Results(i).model;
    subject = Results(i).subject;
    w = Results(i).w;
    dv = Results(i).dv;
    tests = Results(i).tests;
    
    for j = 1:length(ResultsAveraged)
        if strcmp(model, ResultsAveraged(j).model) &&...
                subject == ResultsAveraged(j).subject &&...
                w == ResultsAveraged(j).w &&...
                dv == ResultsAveraged(j).dv
            
            ResultsAveraged(j).tests = ResultsAveraged(j).tests(:,1:9) + tests(:,1:9);
            ResultsAveraged(j).fPos = [ResultsAveraged(j).fPos tests(:,2)];
            ResultsAveraged(j).fSpd = [ResultsAveraged(j).fSpd tests(:,5)];
            ResultsAveraged(j).fAcc = [ResultsAveraged(j).fAcc tests(:,8)];
            ResultsAveraged(j).mPos = [ResultsAveraged(j).mPos tests(:,3)];
            ResultsAveraged(j).mSpd = [ResultsAveraged(j).mSpd tests(:,6)];
            ResultsAveraged(j).mAcc = [ResultsAveraged(j).mAcc tests(:,9)];
            ResultsAveraged(j).n = ResultsAveraged(j).n + 1;
        end
        
        if strcmp(model, ResultsAveraged(j).model) &&...
                ResultsAveraged(j).subject == 0 &&...
                w == ResultsAveraged(j).w &&...
                dv == ResultsAveraged(j).dv
            
            ResultsAveraged(j).tests = ResultsAveraged(j).tests(:,1:9) + tests(:,1:9);
            ResultsAveraged(j).fPos = [ResultsAveraged(j).fPos tests(:,2)];
            ResultsAveraged(j).fSpd = [ResultsAveraged(j).fSpd tests(:,5)];
            ResultsAveraged(j).fAcc = [ResultsAveraged(j).fAcc tests(:,8)];
            ResultsAveraged(j).mPos = [ResultsAveraged(j).mPos tests(:,3)];
            ResultsAveraged(j).mSpd = [ResultsAveraged(j).mSpd tests(:,6)];
            ResultsAveraged(j).mAcc = [ResultsAveraged(j).mAcc tests(:,9)];
            ResultsAveraged(j).n = ResultsAveraged(j).n + 1;
        end
    end  
end

t = (-0.5+1/Hz : 1/Hz : 5)'; % time stamps
% compute confidence interval and divid by n to get the average, n is the
% repetition of each condition
for i = 1:length(ResultsAveraged)
    subject = ResultsAveraged(i).subject;
    n = ResultsAveraged(i).n;
    fPos = ResultsAveraged(i).fPos;
    fSpd = ResultsAveraged(i).fSpd;
    fAcc = ResultsAveraged(i).fAcc;
    mPos = ResultsAveraged(i).mPos;
    mSpd = ResultsAveraged(i).mSpd;
    mAcc = ResultsAveraged(i).mAcc;
    
    if subject == 0
        t_critical = 1.98; % df = n-1, confidence level = 95%, two-tails
    else
        t_critical = 2.262; % df = n-1 = 9, confidence level = 95%, two-tails
    end
    ResultsAveraged(i).data_CI = t_critical/sqrt(n) * [std(fPos,0,2),...
        std(fSpd,0,2),...
        std(fAcc,0,2)];
    ResultsAveraged(i).sim_CI = t_critical/sqrt(n) * [std(mPos,0,2),...
        std(mSpd,0,2),...
        std(mAcc,0,2)];
    ResultsAveraged(i).tests = ResultsAveraged(i).tests ./ n;
    ResultsAveraged(i).tests(:,10) = t;
end


%% Plot average data and simulation by model and subject with confidence interval

iModel = 5; % indicate the model you want to plot by its index in modelList
iSubject = 0; % 0 is the average of all subjects
figure;
iPlot = 1;
for i = 1:length(ResultsAveraged)
    model = ResultsAveraged(i).model;
    subject = ResultsAveraged(i).subject;
    w = ResultsAveraged(i).w;
    dv = ResultsAveraged(i).dv;
    if strcmp(model, modelList(iModel)) && subject == iSubject      
        subplot(3,2,iPlot);
        iPlot = iPlot + 1;
        l = ResultsAveraged(i).tests(:,4);
        f = ResultsAveraged(i).tests(:,5);
        m = ResultsAveraged(i).tests(:,6);
        fCI = ResultsAveraged(i).data_CI(:,2);
        mCI = ResultsAveraged(i).sim_CI(:,2);
        fLow = f - fCI;
        fHigh = f + fCI;
        mLow = m - mCI;
        mHigh = m + mCI;
        
        % confidence intervals
        patch([t;fliplr(t')'], [fLow;fliplr(fHigh')'], 'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'none');
        patch([t;fliplr(t')'], [mLow;fliplr(mHigh')'], 'b', 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none');
        
        hold on;
        plot(t, l, 'k'); % plot the leader
        plot(t, f, 'r'); % plot the follower
        plot(t, m, 'b--'); % plot the model
        title(['w=', num2str(w), '  dv=', num2str(dv), ' ', model{1}]);
        axis([-0.5 6 0.5 1.6]);
        xlabel('Time');
        ylabel('Speed');
%         saveas(gcf,['w=', num2str(w), '  dv=', num2str(dv), ' ', model{1}, '.png']);
    end
end




