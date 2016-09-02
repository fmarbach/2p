% compute and save avg pip rate differences per session
clear;clc
basedir = '/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone/';
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

% read text file with experiment info
fileID = fopen([basedir 'dirList.txt'],'r');
tmp = textscan(fileID,'%s %d');
dirList = tmp{1};
depthList = tmp{2};

thr = 0.5;

for ii = 1:length(dirList)
    [~,sessionName] = fileparts([basedir dirList{ii}]);
    load([basedir dirList{ii} filesep sessionName '.mat'])
    [vec,tr] = dj_trial_extract(dj);
    
    % ensure that we start trial off and end trial off
    tmp = vec.trialOn(find(vec.trialOn==0,1,'first'):find(vec.trialOn==0,1,'last'));
    trialOff = find(diff(tmp)==-1);
    trialOn = find(diff(tmp)==1);

    rF = [];
    rD = [];
    block = [];
    ind = 0;
    for tt = 1:length(trialOn)
        s = trialOn(tt):trialOff(tt);
        if length(s)/500 > thr
            ind = ind+1;
            rF(ind) = sum(vec.pulse(s))/(length(s)/500);
            rD(ind) = sum(vec.ch2_pulse(s))/(length(s)/500);
            block(ind) = mode(vec.ch2_order(s)==12);
        end
    end
    pipDiff(ii) = mean(rF-rD); 
end

save([basedir 'pipDiff.mat'],'pipDiff');

