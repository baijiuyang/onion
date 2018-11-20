clear all;
load Carrot4_data_piped;

% load data_averaged;

%% plot ending speed
figure;

%% plot ending distance
figure;

%% plot leader and follower on each trial
Hz = 90;
for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6*Hz;
    plot(1:length(data_set(i).data(t_start:t_end,3)),data_set(i).data(t_start:t_end,5));
    hold on;
    plot(1:length(data_set(i).data(t_start:t_end,4)),data_set(i).data(t_start:t_end,6));
    d0 = num2str(data_set(i).d0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', dv]);
    axis([0 400 -0.6 0.6]);
%     saveas(fig, [path, 'xRMSEs', '_', num2str(i), '(d0=', d0, ' dv=', dv, ').png']);
    pause(0.2);
    hold off;
    
end

%% plot 1 iteration fitting of 1 model to all data
figure;
trial = 200;
plot(results.test_set(trial).tests(:,4));
hold on;
plot(results.test_set(trial).tests(:,3));

%% plot many trials
figure;
hold on;
for i = 1:length(results.test_set)
    if results.test_set(i).d0 == 8
        plot(results.test_set(i).tests(:,4));
%         plot(results.test_set(i).tests(:,3));
    end
end



%% Pearson's r
Pearson_r = tanh(mean([results.test_set.z_v]));
Pearson_r

%% BIC (wikipedia) form the error of 1 iteration fitting of 1 model to all data
RMSE_v = [Results.RMSE_v];
n = 60; % some sort of degree of freedom. subject*(condition - 1)
for i = 1:length(initParamList)
    k = length(initParamList{i}); % number of parameters
    MSE = mean(RMSE_v(1+(i-1)*706:706+(i-1)*706).^2);
    BIC_wiki(i,1) = n*log(MSE) + k*log(n);
end

%%
for i = 1:length(ResultsSummary)
    ResultsSummary(i).BIC_wiki = BIC_wiki(i);
end


%% BIC (book) form the error of 1 iteration fitting of 1 model to all data
k = length(initParams); % number of parameters
test_length = 361;
n = test_length * n_trials;
RMSE = mean([results.test_set.RMSE_v]);
X = [];
for i = 1:length(results.test_set)
    X = [X; results.test_set(i).tests(:,3)];
end
VAR = var(X);
VARs = zeros(length(results.test_set),1);
for i = 1:length(results.test_set)
    VARs(i) = var(results.test_set(i).tests(:,3));
end
VAR_mean = mean(VARs);

BIC_book = n*log(RMSE) + k*log(n) - (n - k)*log(1 - k/n) + k*log((VAR_mean/RMSE-1)/k);
BIC_book

%% Ploting d, v at the onset of perturbation
distance = [];
velocity = [];
Hz = 90;
for i = 1:length(ExperimentalTrials)
    if ExperimentalTrials(i).dv ~= 0 && ExperimentalTrials(i).dump == 0&&...
            ExperimentalTrials(i).d0 == 1
        manipOnset = int32(ExperimentalTrials(i).manipOnset * Hz);
        curLPos = ExperimentalTrials(i).data(manipOnset,1);
        curFPos = ExperimentalTrials(i).data(manipOnset,2);
        curDist = curLPos - curFPos;
        curFSpd = ExperimentalTrials(i).data(manipOnset,4);
        distance(end+1) = curDist;
        velocity(end+1) = curFSpd;
    end
end
figure;
histogram(distance);
figure;
histogram(velocity);

%% plot speed data against one model
Hz = 90;
x = 1/Hz:1/Hz:5.5;
x = x-0.5;
model = 1;

for i = 1:18
%     (1+(i-1)*706:706+(i-1)*706)
    d0 = num2str(model_ave(i).d0);
    dv = num2str(model_ave(i).dv);
    condition = ['Initial distance ' d0 '(m) leader speed change ' dv '(m/s)'];
    figure;
    hold on;
    fig = plot(x, condition_ave(i).data(:,3), 'k', 'LineWidth',2),
    plot(x, condition_ave(i).data(:,4), 'r', 'LineWidth',2),
    plot(x, model_ave(i).data(:,2), '--r', 'LineWidth',2),
    plot([0 0],[0 2],'--','Color',[0.3 0.3 0.3]),
    text(0, 1.8, '\leftarrow leader speed change'),
    legend('leader','follower',model(5:end)),
    axis([-0.5 5.5 0 2]),
    xlabel('Time(s)'),
    ylabel('Speed(m/s)'),
    title(condition);
%     saveas(fig,['data_vs_', model(5:end), 'p[', num2str(p), ']', '(d0=', d0, ' dv=', dv, ').png']);
end



%% cross-correlation
Hz = 60;
data_set = ExperimentalTrials(1);
for trial = ExperimentalTrials
    if trial.dump == 0 && trial.t_total >= 8.5 && trial.dv ~= 0
        data_set(end+1) = trial;
    end
end
data_set = data_set(2:end);

for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6.0*Hz;
    for j = 1:Hz*3
        j_t_start = t_start + j;
        j_t_end = t_end - j;
        data_set(i).xCorrs(j,1) = corr(data_set(i).data(t_start:j_t_end,3), data_set(i).data(j_t_start:t_end,4));
        data_set(i).xCorrs(j,2) = j;
    end
    [M, I] = max(data_set(i).xCorrs(:,1));
    data_set(i).opt_delay = I;
    data_set(i).opt_xCorr = M;
    data_set(i).xCorr_SD = std(data_set(i).xCorrs(:,1));
    data_set(i).xCorr_range = range(data_set(i).xCorrs(:,1));

end
% save('data_xCorrs.mat', 'data_set');

%% Cross RMSE analysis for finding optimal delay
Hz = 60;
data_set = ExperimentalTrials(1);
for trial = ExperimentalTrials
    if trial.dump == 0 && trial.t_total >= 8.5 && trial.dv ~= 0
        data_set(end+1) = trial;
    end
end
data_set = data_set(2:end);

for i = 1:length(data_set)
    t_start = int32((data_set(i).manipOnset-0.5)*Hz)+1;
    t_end = t_start + 6.0*Hz;
    for j = 1:Hz*5
        j_t_start = t_start + j;
        j_t_end = t_end - j;
        data_set(i).xRMSEs(j,1) = sqrt(mean((data_set(i).data(t_start:j_t_end,3) - data_set(i).data(j_t_start:t_end,4)).^2));
        data_set(i).xRMSEs(j,2) = j;
    end
    [M, I] = min(data_set(i).xRMSEs(:,1));
    data_set(i).opt_delay = I;
    data_set(i).opt_xRMSE = M;
    data_set(i).xRMSE_SD = std(data_set(i).xRMSEs(:,1));
    data_set(i).xRMSE_range = range(data_set(i).xRMSEs(:,1));

end
% save('data_xRMSEs.mat', 'data_set');

%% Histogram of optimal delay
delay = [];
for i = 1:length(data_set)
    if data_set(i).opt_delay < 120 && data_set(i).opt_delay > 3 && data_set(i).d0 == 8 && data_set(i).dv == 0.3
        delay(end+1) = data_set(i).opt_delay/60;
    end
end

% ['Mean = ', num2str(mean(delay))]
% ['Median = ', num2str(median(delay))]
% ['Mode(0.1) = ', num2str(mode(round(delay, 1)))]
edges = -0.05:0.1:2.05;
histogram(delay,edges);
ylim([0 25]);
% histogram(delay,500);
% xticks(0:0.05:3);

%% xRMSE in each trial
path = 'C:\Users\jbai5\OneDrive\First year project\Data analysis\delay\';
for i = 1:length(data_set)
    fig = plot(data_set(i).xRMSEs(:,2),data_set(i).xRMSEs(:,1));   
    d0 = num2str(data_set(i).d0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', dv]);
    axis([0 400 0 0.3]);
%     saveas(fig, [path, 'xRMSEs', '_', num2str(i), '(d0=', d0, ' dv=', dv, ').png']);
    pause(0.2);
    
end

%% cross correlations in each trial
path = 'C:\Users\jbai5\OneDrive\First year project\Data analysis\delay\';
for i = 1:length(data_set)
    fig = plot(data_set(i).xCorrs(:,2),data_set(i).xCorrs(:,1));   
    d0 = num2str(data_set(i).d0);
    dv = num2str(data_set(i).dv);
    title([num2str(i), '  ', d0, ' ', dv]);
    axis([0 400 0 1]);
%     saveas(fig, [path, 'xCorrs', '_', num2str(i), '(d0=', d0, ' dv=', dv, ').png']);
    pause(0.2);
    
end