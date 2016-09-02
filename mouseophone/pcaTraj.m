clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

[~,sessionName]=fileparts(basedir);
load([basedir filesep sessionName '.mat'])
[vec,tr] = dj_trial_extract(dj);

% chosenRois
tmp=dir([basedir filesep  'chosenRois.mat']);
tmp=load([basedir filesep tmp.name]);
S = real(tmp.S);
C = tmp.C;
regmat = dir([basedir filesep  '*00reg.mat']);
regmat = load([basedir filesep regmat.name]);

nextFileVecInd=find(diff(vec.sync)==1);%% find the time at which the sync went up
%%% init stuff 
syncVec=nan(size(vec.sync));
syncFrame=nan(size(S,2),1);
FrameOrigTime=nan(size(S,2),1);

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
    
%     MatchSyncByTime=round(timeNow(1)/10);
%     bhinx=MatchSyncByTime;  %%% Use this one to sunc by pulse timing
    bhinx=ii;    %%% Use this one if no sync pulse are missing
    
    timeBh=vec.t(nextFileVecInd(bhinx-1):nextFileVecInd(bhinx)-1);
    %%% the frame of each behavioral sample in imaging
%     syncVec(nextFileVecInd(ii-1):nextFileVecInd(ii)-1)=interp1(framesNow,timeNow,linspace(framesNow(1),framesNow(end),length(timeBh)));
    
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


%% plot pca trajectories

col = [0.7 0.7 0; 0 0.2 0.8];

%X = C';
X = S'; % time x N
thr = 0.1;
[coeff, score, lat, tsq, expl] = pca(X);

U = coeff(:,1:20); % N x components

%% get Frames of certain trials
c = [2 3 4];
trialOffSamples(1) = 0 * vec.sr; % in samples
trialOffSamples(2) = 0 * vec.sr; % in samples
trialOffset = find(diff(vec.trialOn)==-1);
trialOnset = find(diff(vec.trialOn)==1);
%anchor{2} = find(diff([0; vec.trialOn])==-1 & vec.reward==0);
anchor{1} = find(vec.reward==1 & vec.ch2_order==12);
anchor{2} = find(vec.reward==1 & vec.ch2_order==21);
for ss = 1:length(anchor) % time to align to
    allSamples = {};
    allFrames = {};
    ind = {};
    anchorF = [];
    incr = 0;
    for ii = 1:length(anchor{ss}) % loop over trials
        thisTrial = find(anchor{ss}(ii) > trialOnset,1,'last');
        firstSample = trialOnset(thisTrial);
        if thisTrial < length(trialOnset) % boundary conditions
            nextOnset = trialOnset(thisTrial+1);
        else nextOnset = NaN;
        end
        lastSample = min([nextOnset anchor{ss}(ii)+trialOffSamples(ss) length(vec.trialOn)]);
        trialDur = length(firstSample:anchor{ss}(ii))/vec.sr;
        allSamples{ii} = firstSample:lastSample-1;
        [allFrames{ii},ind{ii}] = unique(round(syncVec(allSamples{ii})));
        
        if ~sum(isnan(allFrames{ii}))>0 && trialDur>thr
            incr = incr+1;
            anchorF = round(syncVec(anchor{ss}(ii)))-allFrames{ii}(1);
            thisX = X(allFrames{ii},:); % time x N
            t = thisX*U; % time x components
            plot3(t(:,c(1)),t(:,c(2)),t(:,c(3)),'LineStyle','none','Marker','o','Color',col(ss,:));
            hold on;
            plot3(t(anchorF,c(1)),t(anchorF,c(2)),t(anchorF,c(3)),'k*');
            plot3(t(1,c(1)),t(1,c(2)),t(1,c(3)),'g*');
            plot3(t(end,c(1)),t(end,c(2)),t(end,c(3)),'r*');
            
%             figure(ss)
%             for pp = 1:20 % loop over components
%                 subplot(4,5,pp);
%                 axis tight
%                 plot(t(:,pp),'Color',col(ss,:));
%                 hold on
%                 plot(t(anchorF(ii),pp),'k*');
% 
%             end
        end
    end

end


%% get same length trial chunks

sortPipRate = 1;
plotsum = 0;
winImg = -23:24; % in frames
inZoneWin = -0.3*vec.sr:0; % 300ms before trial end
inZoneMin = 0.5; % at least fraction in zone in last 300ms
winStatPre = [-round(0.3*imagingSR):0] + abs(winImg(1)); % [sec] window to average for stats
winStatPost = [1:round(0.2*imagingSR)] + abs(winImg(1));

winImgAx=(1000/imagingSR)*winImg'; % x axis in ms
winBeh = round(vec.sr*winImgAx(1)/1000):round(vec.sr*winImgAx(end)/1000);
winBehAx = 1000*winBeh/vec.sr;

order_array=[12 21];
for oo=1:2
    thisInd = vec.reward==1 & vec.ch2_order==order_array(oo);
    rewardFrames{oo} = round(syncVec(thisInd));
    rewardSamples{oo} = find(thisInd);
    badInd = isnan(rewardFrames{oo}) | rewardFrames{oo}<winImg(end);
    rewardFrames{oo}(badInd) = [];
    rewardSamples{oo}(badInd) = [];
    
    % find distractor samples corresponding to rewards in feedback
    rewardSamplesShift{oo} = rewardSamples{oo} + vec.uni_ch2_delay*vec.sr;
    rewardSamplesShift{oo}(rewardSamplesShift{oo} > length(syncVec)) = [];
    rewardFramesShift{oo} = round(syncVec(rewardSamplesShift{oo}));
    badInd = isnan(rewardFramesShift{oo}) | rewardFramesShift{oo}<winImg(end);
    rewardFramesShift{oo}(badInd) = [];
    rewardSamplesShift{oo}(badInd) = [];
    sameBlock = vec.ch2_order(rewardSamplesShift{oo})==order_array(oo);
    switchSamples{oo} = rewardSamplesShift{oo}(~sameBlock);
    switchFrames{oo} = rewardFramesShift{oo}(~sameBlock);
    rewardFramesShift{oo}(~sameBlock) = [];
    rewardSamplesShift{oo}(~sameBlock) = [];
    
    % process unrewarded trials as control
    trialEnd = [0; diff(vec.trialOn)==-1];
    noRewardSamples{oo} = find(trialEnd & ~vec.reward & vec.ch2_order==order_array(oo));
    tmpInd = bsxfun(@plus,noRewardSamples{oo},inZoneWin);
    fracInZone = mean(vec.in_zone(tmpInd),2);
    [fracInZone, inZoneInd] = sort(fracInZone);
    noRewardSamples{oo} = noRewardSamples{oo}(inZoneInd);
    noRewardSamples{oo}(fracInZone < inZoneMin) = [];
    
    noRewardFrames{oo} = round(syncVec(noRewardSamples{oo}));
    badInd = isnan(noRewardFrames{oo}) | noRewardFrames{oo}<winImg(end);
    noRewardFrames{oo}(badInd) = [];
    noRewardSamples{oo}(badInd) = [];
end
rewardSamplesShift{1} = [rewardSamplesShift{1}; switchSamples{2}];
rewardSamplesShift{2} = [rewardSamplesShift{2}; switchSamples{1}];
rewardFramesShift{1} = [rewardFramesShift{1}; switchFrames{2}];
rewardFramesShift{2} = [rewardFramesShift{2}; switchFrames{1}];

% plot pca components
for oo=1:2
    
    % reward responses F1/F2
    rasterR = bsxfun(@plus,rewardFrames{oo},winImg);
    rasterNR = bsxfun(@plus,noRewardFrames{oo},winImg);
    figure(oo)
    for ii = 1:size(rasterInd,1) % loop over trials
        thisX = X(rasterInd(:,ii),:); % time x N
        t = thisX*U; % time x components
        for pp = 1:15 % loop over components
            subplot(3,5,pp);
            axis tight
            plot(t(:,pp),'Color',col(ss,:));
        end
    end
    
    
    % no reward responses F1/F2
    
    raster = thisRoi(rasterInx);
    h4(oo)=subplot(3,4,oo+4+2);
    imagesc(winImgAx,1:size(raster,1),raster)

    
    % distractor responses F1/F2
    rasterInx=bsxfun(@plus,rewardFramesShift{oo},winImg);
    raster=thisRoi(rasterInx);
    h6(oo)=subplot(3,4,2*oo+8);
    imagesc(winImgAx,1:size(raster,1),raster)
    h6(oo).YDir='normal';
    colormap('hot');
    axis tight square;box off
end
Dpre{oo} = mean(raster(:,winStatPre),2);
Dpost{oo} = mean(raster(:,winStatPost),2);





