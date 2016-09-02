clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

[~,sessionName]=fileparts(basedir);
load([basedir filesep sessionName '.mat'])
[vec,tr ] = dj_trial_extract(dj);

foopsiOP = dir([basedir filesep  'chosenRois.mat']); % chosenRois
foopsiOP = load([basedir filesep foopsiOP.name]);
S_df = real(foopsiOP.S);
C_df = foopsiOP.C;
regmat = dir([basedir filesep  '*00reg.mat']);
regmat = load([basedir filesep regmat.name]);

nextFileVecInd=find(diff(vec.sync)==1);% find the time at which the sync went up
% init stuff 
syncVec=nan(size(vec.sync));
syncFrame=nan(size(S_df,2),1);
FrameOrigTime=nan(size(S_df,2),1);

for ii=1:min([size(regmat.results{1},1) length(nextFileVecInd)])
if ii==1 %% first file just init timers
    lastFrame=regmat.results{1}(ii).frameInd(end);
elseif ii==size(regmat.results{1},1)%% if the imaging ran on after behavior stopped (usually not the case)
    framesNow=regmat.results{1}(ii).frameInd+lastFrame;
    timeNow=regmat.results{1}(ii).frameTimeStamps;
    FrameOrigTime(framesNow)=timeNow;
    timeVec(framesNow)=timeNow;
else %% any other one
    framesNow=regmat.results{1}(ii).frameInd+lastFrame;
    timeNow=regmat.results{1}(ii).frameTimeStamps;
    FrameOrigTime(framesNow)=timeNow;
    
    bhinx=ii;    %%% Use this one if no sync pulse are missing
    timeBh=vec.t(nextFileVecInd(bhinx-1):nextFileVecInd(bhinx)-1);
   
    %%% nearest frame number in the clock of the beh computer
    syncVec(nextFileVecInd(bhinx-1):nextFileVecInd(bhinx)-1)=interp1(framesNow,framesNow,linspace(framesNow(1),framesNow(end),length(timeBh)));
    %%% the sample of behavior for each frame
    syncFrame(framesNow)=round(linspace(nextFileVecInd(bhinx-1),nextFileVecInd(bhinx)-1,length(framesNow)));%% time of each frame in bh
    lastFrame=framesNow(end);
    %%% time vector
    timeVec(framesNow)=timeNow;
end
end
imagingSR=inv(median(diff(timeVec)));


% ----------------------------------
%% compute observations 
% y1 (avg trial activity) and y2 (reward period activity)

trialOffSamples = round(0.5 * vec.sr);
trialOffset = find(diff(vec.trialOn)==-1);
trialOnset = find(diff(vec.trialOn)==1);

y1 = zeros(length(trialOffset),size(C_df,1));
y2 = zeros(length(trialOffset),size(C_df,1));
f1_1 = zeros(length(trialOffset),1);
f2_1 = f1_1;
f1_2 = f1_1;
f2_2 = f1_1;
ind = 1;
goodTrials = [];
for ii = 1:length(trialOffset) % loop over trials
    thisTrial = find(trialOffset(ii) > trialOnset,1,'last');
    firstSample1 = trialOnset(thisTrial);
    if thisTrial < length(trialOnset) % boundary condition
        nextOnset = trialOnset(thisTrial+1);
    else
        nextOnset = NaN;
    end
    lastSample2 = min([trialOffset(ii)+trialOffSamples nextOnset length(vec.trialOn)]);
    y1Samp = firstSample1:trialOffset(ii)-1; % trial samples
    y2Samp = trialOffset(ii):lastSample2-1; % reward period samples
    y1Fram = unique(round(syncVec(y1Samp)));
    y2Fram = unique(round(syncVec(y2Samp)));
    y1Fram(isnan(y1Fram)) = [];
    y2Fram(isnan(y2Fram)) = [];
    if ~isempty(y1Fram) && ~isempty(y2Fram)
        y1(ind,:) = mean(S_df(:,y1Fram),2);
        y2(ind,:) = mean(S_df(:,y2Fram),2);
        % get f1 and f2 pip rates
        f1_1(ind) = sum(vec.f1(y1Samp))/(length(y1Samp)/vec.sr);
        f2_1(ind) = sum(vec.f2(y1Samp))/(length(y1Samp)/vec.sr);
        f1_2(ind) = sum(vec.f1(y2Samp))/(length(y2Samp)/vec.sr);
        f2_2(ind) = sum(vec.f2(y2Samp))/(length(y2Samp)/vec.sr);
        goodTrials(ind) = ii;
        ind = ind+1;
    end
end


%% use glmfit to get coefficients b
% y is n x 1; observations, y1==trial activity, y2==reward period activity
% X is n x p (n trials, p predictors)
% here we include trials of one session (max available for one cell)
% predictors: 
%   1) avg pip rate f1
%   2) avg pip rate f2
%   3) reward 
%   4) block 
%   5) transfer fct 
%   6) trial dur
x1_1 = f1_1(goodTrials);
x2_1 = f2_1(goodTrials);
x1_2 = f1_2(goodTrials);
x2_2 = f2_2(goodTrials);
x3 = tr.outcome(goodTrials);
x4 = tr.ch2Order(goodTrials)==12;
x5 = tr.level(goodTrials);
x6 = tr.rt(goodTrials);

X1 = [x1_1 x2_1 x3 x4 x5 x6];
X2 = [x1_2 x2_2 x3 x4 x5 x6];

b1 = zeros(size(C_df,1),size(X1,2)+1);
b2 = zeros(size(C_df,1),size(X2,2)+1);
for rr = 1:size(C_df,1)
    b1(rr,:) = glmfit(X1,y1(goodTrials,rr),'normal');%,'link','logit');
    b2(rr,:) = glmfit(X2,y2(goodTrials,rr),'normal');
end






