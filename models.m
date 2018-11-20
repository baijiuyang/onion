function [mPos, mSpd, mAcc] = models(model, p, delay, lPos, lSpd, pStart, vStart, inputHz, outputHz)
if strcmp(model, 'nullModel')
    [mPos, mSpd, mAcc] = nullModel(pStart, vStart, lSpd, inputHz);
    
elseif strcmp(model,'speedModel')
    c = p(1);
    [mPos, mSpd, mAcc] = speedModel(pStart, vStart, inputHz,outputHz, lSpd, c, delay);
    
elseif strcmp(model, 'distanceModel')
    c = p(1);
    [mPos, mSpd, mAcc] = distanceModel(pStart, vStart, inputHz,outputHz, lPos, c, delay);
    
elseif strcmp(model, 'expansionModel')
    b = p(1);
    w = 0.2; % the width of the target pole is 0.6 meter
    [mPos, mSpd, mAcc] = expansionModel(pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay);
    
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
    
elseif strcmp(model, 'expansion2Model')
    b = p(1);
    w = 0.6; % the width of the target pole is 0.6 meter
    [mPos, mSpd, mAcc] = expansion2Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, delay);
    
elseif strcmp(model, 'expansion3Model')
    b = p(1);
    L = p(2);
    w = 1.2; % the width of the target pole is 0.6 meter
    [mPos, mSpd, mAcc] = expansion3Model( pStart, vStart,inputHz,outputHz, lPos, lSpd, b, w, L, delay);
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
% alpha.: the changing rate of visual angle
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
% alpha.: the changing rate of visual angle
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
% alpha.: the changing rate of visual angle
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


end