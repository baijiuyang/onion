clear all;
load Onion_data_piped;

%% plot all data
figure;
hold on;
Hz = 90;
for j = 1:length(ExperimentalTrials)
    if ExperimentalTrials(j).dump ~= 1 %&& ExperimentalTrials(j).t_total > 8.5 && ...
%             ExperimentalTrials(j).dv ~= 0
        manipOnset = ExperimentalTrials(j).manipOnset;
        tStart = int32((manipOnset-0.5)*Hz) + 1;
        tEnd = int32((manipOnset + 5)*Hz);
        x = ExperimentalTrials(j).data(tStart:tEnd,7)-manipOnset;
        y = ExperimentalTrials(j).data(tStart:tEnd,4);
%         if sum(y >= 1.9)>0
%             j
%         end
        plot(x,y);
        axis([-2 5 0 2]);
    end
end


%% plot filtered freewalk 
for j = 1:length(freewalk)
    hold on;
    plot(freewalk(j).data(:,4), freewalk(j).data(:,2));
end

%% plot input traj
figure;
hold on;
for j = 1:length(ExperimentalTrials)
    if ExperimentalTrials(j).subject >= 1 && ExperimentalTrials(j).dump ~= 1  && ExperimentalTrials(j).dv ~= 0
        x = ExperimentalTrials(j).data(:,7);
        y = ExperimentalTrials(j).data(:,3);
        plot(x,y);
        axis([0 12 0 2]);
    end
end


%% plot rotated and interpolated but unfiltered data

figure;
hold on;
for j = 1:length(ExperimentalTrials)  
    if ExperimentalTrials(j).dump ~= 1 && ExperimentalTrials(j).subject >= 1 
        speed = diff(ExperimentalTrials(j).inter_traj(:,4));
        speed = [speed; speed(end)];
        speed = speed * 90;
        plot(ExperimentalTrials(j).inter_traj(:,5), speed);
%         axis([0 12 -5 5]);
    end
end

%% plot filtered data by conditions in different graphs
figure;
for j = 1:length(ExperimentalTrials)
    subject = ExperimentalTrials(j).subject;
    x = ExperimentalTrials(j).data(:,7);
    y = ExperimentalTrials(j).data(:,4);
    dump = ExperimentalTrials(j).dump;
    w = ExperimentalTrials(j).w;
    dv = ExperimentalTrials(j).dv;
    if dump ~= 1 && subject >= 1
        if dv > 0 && w == 0.2
            subplot(2,3,1);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=0.3 w=0.2');
            axis([0 12 0 2]);
        end
        if dv > 0 && w == 0.6
            subplot(2,3,2);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=0.3 w=0.6');
            axis([0 12 0 2]);
        end
        if dv > 0 && w == 1
            subplot(2,3,3);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=0.3 w=1');
            axis([0 12 0 2]);
        end
        if dv < 0 && w == 0.2
            subplot(2,3,4);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=-0.3 w=0.2');
            axis([0 12 0 2]);
        end
        if dv < 0 && w == 0.6
            subplot(2,3,5);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=-0.3 w=0.6');
            axis([0 12 0 2]);
        end
        if dv < 0 && w == 1
            subplot(2,3,6);
            hold on;
            plot(x,y);
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title('dv=-0.3 w=1');
            axis([0 12 0 2]);
        end
    end
end


%% plot for all subject in individual graph
figure;
for i = 1:12
    subplot(3,4,i);
    hold on;
    for j = 1:length(ExperimentalTrials) 
        if ExperimentalTrials(j).dump ~= 1 && ExperimentalTrials(j).subject == i %&& ExperimentalTrials(j).dv ~= 0
            plot(ExperimentalTrials(j).data(:,7), ExperimentalTrials(j).data(:,4));
            xlabel('Time(s)');
            ylabel('Speed(m/s)');
            title(['subject ' num2str(i)]);
            axis([0 12 0 2]);
        end
    end
end


%% histogram of trial length
t = [];
Hz = 90;
figure;
hold on;
for j = 1:length(ExperimentalTrials)
    if ExperimentalTrials(j).dump ~= 1 %ExperimentalTrials(i).w == 0.2 && ExperimentalTrials(i).dv == -0.3
        t(end+1) = ExperimentalTrials(j).t_total - ExperimentalTrials(j).manipOnset;
    end
end
histogram(t,10);
axis([0 13 0 500]);

%% trial info check (trial length analysis)
count = 0;
for j = 1:length(ExperimentalTrials)
    subject = ExperimentalTrials(j).subject;
    trial = ExperimentalTrials(j).trial;
    w = ExperimentalTrials(j).w;
    dv = ExperimentalTrials(j).dv;
    dump = ExperimentalTrials(j).dump;
    data = ExperimentalTrials(j).data;
    t = ExperimentalTrials(j).t_total;
    manipOnset = ExperimentalTrials(j).manipOnset;
    
    if dump ~= 1 && t - manipOnset > 5.5 && dv ~= 0
        count = count + 1;
        [w dv]
    end
end
count


%% plot orientation

% load data_orientation;
figure;
hold on;
for j = 1:length(orientation)
    % plot yaw
    y = orientation(j).orientation(:,1);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('yaw');
    ax = gca;
    ax.FontSize = 20;
end

figure;
hold on;
for j = 1:length(orientation)
    % plot yaw
    y = orientation(j).orientation(:,2);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('pitch');
    ax = gca;
    ax.FontSize = 20;  
end
figure;
hold on;
for j = 1:length(orientation)
    % plot yaw
    y = orientation(j).orientation(:,3);
    x = 1:length(y);
    x = x';
    plot(x,y);
    title('roll');
    ax = gca;
    ax.FontSize = 20;
end

%% Check frame rate
FR = [];
for j = 1:length(ExperimentalTrials)
    FR = [FR;(ExperimentalTrials(j).raw_traj(2:end,5) - ExperimentalTrials(j).raw_traj(1:end-1,5)).^(-1)];
end

histogram(FR);
xticks(0:5:100);
xlabel('Hz');
ylabel('Number of frames');
axis([20 100 0 200000]),...
ax = gca;
ax.FontSize = 20;

