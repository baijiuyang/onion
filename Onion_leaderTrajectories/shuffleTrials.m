function [shuffledNfo,shuffledPos,shuffledSpd,shuffledHdn] = shuffleTrials(Nfo, Pos, Spd, Hdn)
% This function will take four array and produce four shuffled array

nTrial = length(Nfo);
shuffledNfo = cell(nTrial,1);
shuffledPos = cell(nTrial,1);
shuffledSpd = cell(nTrial,1);
shuffledHdn = cell(nTrial,1);

randomTrial = randperm(nTrial);

for iTrial = 1:nTrial
    shuffledNfo{iTrial} = Nfo{randomTrial(iTrial)};
    shuffledPos{iTrial} = Pos{randomTrial(iTrial)};
    shuffledSpd{iTrial} = Spd{randomTrial(iTrial)};
    shuffledHdn{iTrial} = Hdn{randomTrial(iTrial)};
end