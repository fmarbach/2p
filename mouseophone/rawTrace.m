clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

[~,sessionName]=fileparts(basedir);

load([basedir filesep sessionName '.mat'])

[vec,tr ]=dj_trial_extract(dj);

% chosenRois
foopsiOP=dir([basedir filesep  'chosenRois.mat']);
foopsiOP=load([basedir filesep foopsiOP.name]);
S_df=real(foopsiOP.S);
C_df=foopsiOP.C;
% foopsiOP=dir([basedir filesep  'ds*']);
% b=load([basedir filesep foopsiOP.name]);

regmat=dir([basedir filesep  '*00reg.mat']);
regmat=load([basedir filesep regmat.name]);

nextFileVecInd=find(diff(vec.sync)==1);%% find the time at which the sync went up
%%% init stuff 
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


% ----------------------------------
%% plot raw data with blocks over it

trialOffSamples(1) = 0 * vec.sr; % in samples
trialOffSamples(2) = 0.5 * vec.sr; % in samples
trialOffSamples(3) = 0 * vec.sr; % in samples
trialOffSamples(4) = 0.5 * vec.sr; % in samples
trialOffset = find(diff(vec.trialOn)==-1);
trialOnset = find(diff(vec.trialOn)==1);
anchor{1} = trialOffset;
anchor{2} = trialOffset;
anchor{3} = find(vec.reward==1);
anchor{4} = find(vec.reward==1);
titles = {'all trials; trial only','all trials; +R period',...
          'R trials; trial only','R trials; +R period',...
          'trial off','full trace'};
for ss = 1:length(anchor)+2
    allSamples = {};
    lastSample = [];
    allFrames = {};
    ind = {};
    if ss <= length(anchor)
        for ii = 1:length(anchor{ss})
            thisTrial = find(anchor{ss}(ii) > trialOnset,1,'last');
            firstSample = trialOnset(thisTrial);
            if thisTrial < length(trialOnset) % boundary conditions
                nextOnset = trialOnset(thisTrial+1);
            else
                nextOnset = NaN;
            end
            lastSample = min([nextOnset anchor{ss}(ii)+trialOffSamples(ss) length(vec.trialOn)]);
            allSamples{ii} = firstSample:lastSample-1;
        end
        tmp = cell2mat(allSamples);
    elseif ss == length(anchor)+1
        tmp = find(vec.trialOn==0);
    elseif ss == length(anchor)+2
        tmp = 1:length(vec.trialOn);
    end
    
    [allFrames,ind] = unique(round(syncVec(tmp)));
    blockSamples = vec.ch2_order(tmp)==12;
    blockSamples = blockSamples(ind);
    blockSamples(isnan(allFrames)) = [];
    allFrames(isnan(allFrames)) = [];    
    
    % prepare block shading
    x1 = find(diff(blockSamples)==1);
    x2 = find(diff(blockSamples)==-1);
    if x1(1) > x2(1) % started with block==12
        x1 = [1; x1];
    end
    if x1(end) > x2(end) % ended with block==12
        x2 = [x2; length(allFrames)];
    end
    xvec = [1:length(allFrames)] / imagingSR; % xaxis in [s]
    
    out(ss).xvec = xvec;
    out(ss).x1 = x1;
    out(ss).x2 = x2;
    out(ss).block = blockSamples;
    out(ss).frames = allFrames;
end

% plot things
for rr = 1:size(C_df,1)
    f1=figure(1);clf;
    f1.Units = 'inches';
    f1.Position = [2 2 8 15];
    f1.Name = ['roi ' num2str(rr)];
    
    for ss = 1:length(out)
        h(ss) = subplot(length(out),1,ss);
        thisRoi = S_df(rr,out(ss).frames);
        xvec = out(ss).xvec;
        x1 = out(ss).x1;
        x2 = out(ss).x2;
        thisMax = max(thisRoi);
        thisMin = min(thisRoi);

        for ii = 1:length(out(ss).x1)
            rectangle('Position',[xvec(x1(ii)) thisMin xvec(x2(ii)-x1(ii)) thisMax-thisMin],...
                      'EdgeColor','none','FaceColor',0.85*[1 1 1])
            hold on
        end
        plot(xvec,thisRoi);
        title(titles{ss});
        axis tight;
        
    end
    eqAx('y',h(1:4));
%     formatFigure(hf, ha, textSize, lineWidth, markerSize, xlab, ylab, ...
%                            xlimit, ylimit, tickDirOut, pbasp, figureSize, boxON, whiteBckg)
    fname = [basedir filesep 'plots\roi' num2str(rr) '.jpg'];
    saveas(f1,fname);
    %pause
    
end

    