%% new analysis with foopsi

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

% save trial durations, outcome and %success
fname = [basedir filesep 'trial.mat'];
tDur = tr.rt;
outcome = tr.outcome;
pSuccess = sum(outcome)/length(outcome);
save(fname,'tDur','pSuccess','outcome');


%% simple psth around reward times
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

limTrialNum=min(cellfun(@length,rewardFrames));
% limTrialNum=80;

for rr=1:size(C_df,1)
    if plotsum
    f1=figure(1);clf;
    figureSet(f1,6,6);
    end
    for oo=1:2
        roinow=S_df(rr,:);
        
        % pips raster F1/F2 reward
        rasterInx = bsxfun(@plus,rewardSamples{oo},winBeh);
        raster = vec.(['f' num2str(oo)])(rasterInx) - 1;
        raster(raster==-1) = NaN; % for raster plot appearance
        rasterStepped=bsxfun(@plus,raster,1*[1:size(raster,1)]')';
        if plotsum
        h1(oo)=subplot(3,4,oo);
        plot(winBehAx,rasterStepped,'k.')
        axis tight square;box off
        title(['pips F' num2str(oo) '-R']);
        end
        
        % pips raster F1/F2 no reward
        rasterInx = bsxfun(@plus,noRewardSamples{oo},winBeh);
        raster = vec.(['f' num2str(oo)])(rasterInx) - 1;
        raster(raster==-1) = NaN; % for raster plot appearance
        rasterStepped=bsxfun(@plus,raster,1*[1:size(raster,1)]')';
        if plotsum
        h2(oo)=subplot(3,4,oo+2);
        plot(winBehAx,rasterStepped,'k.')
        axis tight square;box off
        title(['pips F' num2str(oo) '-noR']);
        end
        
        % reward responses F1/F2
        rasterInx=bsxfun(@plus,rewardFrames{oo},winImg);
        raster=roinow(rasterInx);
        if plotsum
        h3(oo)=subplot(3,4,oo+4);
        imagesc(winImgAx,1:size(raster,1),raster)
        h3(oo).YDir='normal';
        colormap('hot');
        axis tight square;box off
        title(['F' num2str(oo) '-R']);
        end
        Fpre{oo} = mean(raster(:,winStatPre),2);
        Fpost{oo} = mean(raster(:,winStatPost),2);
        
        % no reward responses F1/F2
        rasterInx = bsxfun(@plus,noRewardFrames{oo},winImg);
        raster = roinow(rasterInx);
        if plotsum
        h4(oo)=subplot(3,4,oo+4+2);
        imagesc(winImgAx,1:size(raster,1),raster)
        h4(oo).YDir='normal';
        colormap('hot');
        axis tight square;box off
        title(['F' num2str(oo) '-noR']);
        end
        
        % distractor pips for shifted reward
        rasterInx = bsxfun(@plus,rewardSamplesShift{oo},winBeh);
        raster = vec.ch2_pulse(rasterInx) - 1;
        raster(raster==-1) = NaN; % for raster plot appearance
        rasterStepped=bsxfun(@plus,raster,1*[1:size(raster,1)]')';
        if plotsum
        h5(oo)=subplot(3,4,2*oo-1+8);
        plot(winBehAx,rasterStepped,'k.')
        axis tight square;box off
        if oo==1 title('D pips f2');
        else title('D pips f1'); end
        end
        
        % distractor responses F1/F2
        rasterInx=bsxfun(@plus,rewardFramesShift{oo},winImg);
        raster=roinow(rasterInx);
        if plotsum
        h6(oo)=subplot(3,4,2*oo+8);
        imagesc(winImgAx,1:size(raster,1),raster)
        h6(oo).YDir='normal';
        colormap('hot');
        axis tight square;box off
        if oo==1 title('D f2');
        else title('D f1'); end
        end
        Dpre{oo} = mean(raster(:,winStatPre),2);
        Dpost{oo} = mean(raster(:,winStatPost),2);
        
        % traces
%         rasterInx=bsxfun(@plus,rewardFrames{oo}(1:limTrialNum),winImg);
%         raster=roinow(rasterInx);
%         rasterStepped=bsxfun(@plus,raster,0.15*[1:size(raster,1)]')';
%         h2(oo)=subplot(2,2,oo);
%         plot(winImgAx,rasterStepped,'k');
        %colormap(flipud(colormap('gray')))
        
        % confidence intervals
        Dpre_ci{oo} = bootci(1000,@mean,Dpre{oo});
        Dpost_ci{oo} = bootci(1000,@mean,Dpost{oo});
        Fpre_ci{oo} = bootci(1000,@mean,Fpre{oo});
        Fpost_ci{oo} = bootci(1000,@mean,Fpost{oo});
    end

    % -- plotting --
    if plotsum
    eqAxC([h3 h4 h6]);
    eqAx('both',h2);
    figure(2);
    % pre reward
    h7(1) = subplot(1,2,1);
    y = mean(Fpre{1}); % feedback pips freq1
    errorbar(1,y,abs(Fpre_ci{1}(1)-y),abs(Fpre_ci{1}(2)-y),'ko'); % first lower, then upper ci
    hold on
    %plot(2*ones(1,length(Fpre{1})),Fpre{1},'ko')
    y = mean(Dpre{2}); % distractor pips freq1
    errorbar(2,y,abs(Dpre_ci{2}(1)-y),abs(Dpre_ci{2}(2)-y),'ro'); 
    y = mean(Fpre{2}); % feedback pips freq2
    errorbar(3,y,abs(Fpre_ci{2}(1)-y),abs(Fpre_ci{2}(2)-y),'ko'); 
    y = mean(Dpre{1}); % distractor pips freq2
    errorbar(4,y,abs(Dpre_ci{1}(1)-y),abs(Dpre_ci{1}(2)-y),'ro'); 
    hold off
    % post reward
    h7(2) = subplot(1,2,2);
    y = mean(Fpost{1}); % feeback pips freq1
    errorbar(1,y,abs(Fpost_ci{1}(1)-y),abs(Fpost_ci{1}(2)-y),'ko'); 
    hold on
    y = mean(Dpost{2}); % distractor pips freq1
    errorbar(2,y,abs(Dpost_ci{2}(1)-y),abs(Dpost_ci{2}(2)-y),'ro'); 
    y = mean(Fpost{2}); % feeback pips freq1
    errorbar(3,y,abs(Fpost_ci{2}(1)-y),abs(Fpost_ci{2}(2)-y),'ko'); 
    y = mean(Dpost{1}); % distractor pips freq1
    errorbar(4,y,abs(Dpost_ci{1}(1)-y),abs(Dpost_ci{1}(2)-y),'ro');
    hold off
    eqAx('both',h7);
    end
    
    % do stats
    psum(rr,1) = ranksum(Fpre{1},Fpre{2}); % F1 vs F2
    psum(rr,2) = ranksum(Dpre{1},Dpre{2}); % D1 vs D2
    psum(rr,3) = ranksum(Fpre{1},Dpre{2}); % freq1 F vs D
    psum(rr,4) = ranksum(Fpre{2},Dpre{1}); % freq2 F vs D
        
    if plotsum
%     saveas(f1,['figureTemp' filesep '4RastersBar_' num2str(rr) '.pdf'])
    formatFigure(f1, h1, 10, 1.5, 0, 0, 0, ...
                           0, 0, 1, 0, [1 1 8 6], 1, 1)
    %eqAx('both',hs);
    fname = [basedir filesep 'plots\raster_roi' num2str(rr)];
    %pause
    printpdf(f1,fname);
    psum(rr,:) < 0.05/(size(C_df,1)*2)
    %pause
    end
end


%% activity in function of pips/S, feedback / distractor

binSize = 100; % bin in samples (20 samples = 40ms)
ploton = 1;
%pipSecBinVec = [0 24 44 64 84 inf];
%pipSecBinCenters = [12 34 54 74 92];
pipSecBinVec = [-inf 25 50 75 inf];
pipSecBinCenters = [12.5 37.5 62.5 87.5];
ind.F1 = vec.trialOn & vec.ch2_order==12 & ~isnan(syncVec);
ind.D1 = vec.trialOn & vec.ch2_order==21 & ~isnan(syncVec);
ind.F2 = vec.trialOn & vec.ch2_order==21 & ~isnan(syncVec);
ind.D2 = vec.trialOn & vec.ch2_order==12 & ~isnan(syncVec);
pips.F1 = vec.f1(ind.F1);
pips.D1 = vec.f1(ind.D1);
pips.F2 = vec.f2(ind.F2);
pips.D2 = vec.f2(ind.D2);

fields = fieldnames(pips);
for rr = 1:size(C_df,1)
    if ploton
    f2=figure(2);clf; 
    figureSet(f2,4,8);
    f1=figure(1);clf; 
    figureSet(f1,4,8);
    end
    for ii = 1:length(fieldnames(pips))
        thisPips = pips.(fields{ii});
        binVec = 1:binSize:length(thisPips);
        thisPips = thisPips(1:binVec(end)-1); % omit samples after last bin
        [~,idx] = histc(1:length(thisPips),binVec);
        nPips  = accumarray(idx',thisPips); % number of pips / bin
        frameInd = round(syncVec(ind.(fields{ii})));
        frameInd = frameInd(1:binVec(end)-1); % omit last incomplete bin
        [frameInd, IA] = unique(frameInd); % there are many beh samples for one imaging frame
        idx_f = idx(IA); % 1=all frames belonging to bin1, 2=frames belonging to bin2, etc...
        thisFrames = S_df(rr,frameInd); % all frames for current condition

        nEvents = accumarray(idx_f',thisFrames,[],@mean); % avg events / bin
        
        % bin the pip/s bin into X bins
        nPipsSec = nPips / (binSize/vec.sr); % convert to pips/s
        [~,pipbins] = histc(nPipsSec,pipSecBinVec); % assign bins given by pipSecBinVec
        
        nEventsCell{ii} = accumarray(pipbins,nEvents,[],@(x) {x}); % given cell contains all time bins that have a certain pips/s
        nf_mean{ii} = cellfun(@mean,nEventsCell{ii});
        for bb = 1:length(nEventsCell{ii}) % bootci does not take cell arrays
            nf_ci{ii}(:,bb) = bootci(1000,@mean,nEventsCell{ii}{bb});
        end
        
        %nf_std{ii} = accumarray(np+1,nf,[],@std);
        xvec{ii} = pipSecBinCenters(unique(pipbins));
        
        % stats: is this roi/condition modulated by pip rate?
        p1(rr,ii) = kruskalwallis(nEvents,pipbins,'off');   
        
        if ploton
        figure(f1); hs(ii) = subplot(1,4,ii);
        scatter(pipbins,nEvents);
        end
    end

    % stats: ranksum test between matching pipSec bins feedback vs. distractor
    for pp = 1:length(pipSecBinCenters)
        p2(1,rr,pp) = ranksum(nEventsCell{1}{pp},nEventsCell{2}{pp}); % freq 1
        p2(2,rr,pp) = ranksum(nEventsCell{3}{pp},nEventsCell{4}{pp}); % freq 2
    end
    % ranksum between all pip bins feedback vs. distractor
    p3(rr,1) = ranksum(cell2mat(nEventsCell{1}),cell2mat(nEventsCell{2}));
    p3(rr,2) = ranksum(cell2mat(nEventsCell{3}),cell2mat(nEventsCell{4}));
    
    % stats: set significance level with multiple comparisons
    sigThr_audOnly = 0.01 / (2*length(pipSecBinCenters)); % account for multiple comparisons
    sigThr_pipMod = 0.01 / length(fieldnames(pips));
    if ploton
    sigPipMod = p1(rr,:) < sigThr_pipMod
    sigDistrF1 = squeeze(p2(1,rr,:) < sigThr_audOnly)'
    sigDistrF2 = squeeze(p2(2,rr,:) < sigThr_audOnly)'
    psumDistr = [psum(rr,3) psum(rr,4)] < (0.01/2)
    sigPipMod1 = sigPipMod(1:2);
    sigPipMod2 = sigPipMod(3:4);
    xvec1 = xvec{1}(1:2);
    xvec2 = xvec1;
    
    % plot for each roi
    figure(f2);
    colorF = [0.7 0.6 0.1];
    colorD = [0.4 0.3 0.3];
    h1 = subplot(121);
    y = nf_mean{1}';
    shadedErrorBar(xvec{1},y,[abs(nf_ci{1}(2,:)-y);abs(nf_ci{1}(1,:)-y)],{'Color',colorF},1);
    hold on
    y = nf_mean{2}';
    shadedErrorBar(xvec{2},y,[abs(nf_ci{2}(2,:)-y);abs(nf_ci{2}(1,:)-y)],{'Color',colorD},1);
    xvec1(sigPipMod1==0) = [];
    sigPipMod1(sigPipMod1==0) = [];
    xvec{2}(sigDistrF1==0) = [];
    sigDistrF1(sigDistrF1==0) = [];
    y = min([nf_ci{1}(1,:) nf_ci{2}(1,:)]);
    plot(xvec1+10,y*sigPipMod1,'k*');
    plot(xvec{2},y*sigDistrF1,'r*');
    hold off;
    title(['roi' num2str(rr) ' freq1']);
    
    h2 = subplot(122);
    y = nf_mean{3}';
    shadedErrorBar(xvec{3},y,[abs(nf_ci{3}(2,:)-y);abs(nf_ci{3}(1,:)-y)],{'Color',colorF},1);
    hold on
    y = nf_mean{4}';
    shadedErrorBar(xvec{4},y,[abs(nf_ci{4}(2,:)-y);abs(nf_ci{4}(1,:)-y)],{'Color',colorD},1);
    xvec2(sigPipMod2==0) = [];
    sigPipMod2(sigPipMod2==0) = [];
    xvec{3}(sigDistrF2==0) = [];
    sigDistrF2(sigDistrF2==0) = [];
    y = min([nf_ci{3}(1,:) nf_ci{4}(1,:)]);
    plot(xvec2+10,y*sigPipMod2,'k*');
    plot(xvec{3},y*sigDistrF2,'r*');
    hold off;
    title(['roi' num2str(rr) ' freq2']);
    
    eqAx('both',[h1 h2]);
%     formatFigure(f2, [h1 h2], textSize, lineWidth, markerSize, xlab, ylab, ...
%                            xlimit, ylimit, tickDirOut, pbasp, figureSize, boxON, whiteBckg)
    formatFigure(f2, [h1 h2], 12, 1.5, 0, 'pips/S', 'calcium events [a.u.]', ...
                           0, 0, 1, 0, [1 1 6 4], 1, 1)
    %eqAx('both',hs);
    fname = [basedir filesep 'plots\pipMod_roi' num2str(rr)];
    %pause
    printpdf(f2,fname);
    end
end

% save p values
sigThr_sum = 0.01/2;
psumH = psum < sigThr_sum;
p1H = p1 < sigThr_pipMod;
p2H = p2 < sigThr_audOnly;
fname = [basedir filesep 'pvalues.mat'];
save(fname,'psum','psumH','p1H','p1','p2','p2H','sigThr_pipMod','sigThr_audOnly','sigThr_sum');

% plot pie chart
figure;
p2H_F1 = sum(squeeze(p2H(1,:,:)),2);
p2H_F2 = sum(squeeze(p2H(2,:,:)),2);
aud = sum(p1H,2) & ~p2H_F1 & ~p2H_F2;
feed = p2H_F1 | p2H_F2;
nomod = ~aud & ~feed;
pie([sum(aud) sum(feed) sum(nomod)],{'aud','feed','nomod'});

figure;
feed = psumH(:,3) | psumH(:,4);
nomod = ~feed;
pie([sum(feed) sum(nomod)],{'feed','nomod'});

figure;
p3H = p3 < 0.01/2;
feed = p3H(:,1) | p3H(:,2);
nomod = ~feed;
pie([sum(feed) sum(nomod)],{'feed','nomod'});
