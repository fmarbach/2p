
clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));
[~,sessionName]=fileparts(basedir);
load([basedir filesep sessionName '.mat'])
[vec,tr] = dj_trial_extract(dj);

% chosenRois
foopsiOP=dir([basedir filesep  'chosenRois.mat']);
foopsiOP=load([basedir filesep foopsiOP.name]);
S_df=real(foopsiOP.S);
C_df=foopsiOP.C;
regmat=dir([basedir filesep  '*00reg.mat']);
regmat=load([basedir filesep regmat.name]);
nextFileVecInd=find(diff(vec.sync)==1);%% find the time at which the sync went up
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
    bhinx=ii;
    timeBh=vec.t(nextFileVecInd(bhinx-1):nextFileVecInd(bhinx)-1);
    syncVec(nextFileVecInd(bhinx-1):nextFileVecInd(bhinx)-1)=interp1(framesNow,framesNow,linspace(framesNow(1),framesNow(end),length(timeBh)));
    syncFrame(framesNow)=round(linspace(nextFileVecInd(bhinx-1),nextFileVecInd(bhinx)-1,length(framesNow)));%% time of each frame in bh
    lastFrame=framesNow(end);
    timeVec(framesNow)=timeNow;
end
end
imagingSR=inv(median(diff(timeVec)));


% =====================================================
% bin trial off data, classify into low rate / high rate
% avoid 500ms just after trial
% avoid short ITIs
% then use sliding window to compute pip/s vs. activity

w = 0.2; % sliding window in [s]
binSize = 5; % in [s]
%pipSecBinVec = [0 24 44 64 84 inf];
%pipSecBinCenters = [12 34 54 74 92];
pipSecBinVec = [-inf 25 50 75 inf];
pipSecBinCenters = [12.5 37.5 62.5 87.5];

% find ITI samples
off_offset = 0.5; % in [s]
on_offset = 0.3;
tmp = vec.trialOn(find(vec.trialOn==0,1,'first'):find(vec.trialOn==0,1,'last'));
trialOff = find(diff(tmp)==-1); % ensure that we start trial off and end trial off
trialOn = find(diff(tmp)==1);
trialOff(end) = []; % we are interested in the ITI, so Off-->On
trialOn(1) = [];
ITIsamples = ~vec.trialOn;
for ii = 1:length(trialOff) % extend trialOn period on both ends to avoid reward/trial start
    ITIsamples(trialOff:trialOff+off_offset) = 0;
    ITIsamples(trialOn:trialOn+on_offset) = 0;
end
ind.D1 = ITIsamples & vec.ch2_order==21 & ~isnan(syncVec);
ind.D2 = ITIsamples & vec.ch2_order==12 & ~isnan(syncVec);
pips.D1 = vec.f1(ind.D1);
pips.D2 = vec.f2(ind.D2);

fields = fieldnames(pips);
f2 = figure(2);
for rr = 1:size(C_df,1)
    cc = 1;
    for ii = 1:length(fieldnames(pips)) % loop over types of ITI
        thisPips = pips.(fields{ii});
        thisFrames = round(syncVec(ind.(fields{ii})));
        binVec = 1:binSize*vec.sr:length(thisPips);
        thisPips = thisPips(1:binVec(end)-1); % omit samples after last bin
        [~,idx] = histc(1:length(thisPips),binVec);
        pipRate  = accumarray(idx',thisPips); % number of pips / bin
        pipStd  = accumarray(idx',thisPips,[],@(x) std(diff(find(x)))); % std of IPIs / bin
        lowStd = pipStd<nanmean(pipStd);
        lowRate = pipRate<mean(pipRate) & pipRate>0;
        highRate = pipRate>mean(pipRate);

        lR = [];
        hR = [];
        p(1).ind = [];
        p(2).ind = [];
        for jj = 1:length(binVec)-1 % concatenate similar bins
            if lowStd(jj) && lowRate(jj)
                p(1).ind = [p(1).ind binVec(jj):binVec(jj+1)-1];
            elseif lowStd(jj) && highRate(jj)
                p(2).ind = [p(2).ind binVec(jj):binVec(jj+1)-1];
            end
        end

        % get corresponding imaging frames indeces
        frameInd = round(syncVec(ind.(fields{ii})));
        frameInd(binVec(end):end) = []; % get rid of last incomplete bin

        % pass sliding window over those samples and calculate 
        % mean activity vs. pipRate
        for pp = 1:length(p) % loop over high/low rate epochs
            relevantPips = thisPips(p(pp).ind);
            relevantFrameInd = thisFrames(p(pp).ind);
            binVec = 1:w*vec.sr:length(relevantPips);
            relevantPips = relevantPips(1:binVec(end)-1); % omit samples after last bin
            [~,idx] = histc(1:length(relevantPips),binVec);
            nPips  = accumarray(idx',relevantPips) / w; % number of pips / bin
            
            relevantFrameInd = relevantFrameInd(1:binVec(end)-1); % omit last incomplete bin
            [relevantFrameInd, IA] = unique(relevantFrameInd); % there are many beh samples for one imaging frame
            idx_f = idx(IA); % 1=all frames belonging to bin1, 2=frames belonging to bin2, etc...
            relevantFrames = S_df(rr,relevantFrameInd); % all frames for current condition
            nEvents = accumarray(idx_f',relevantFrames,[],@mean); % avg events / bin

            % bin the pip/s bins into X bins
            [~,pipbins] = histc(nPips,pipSecBinVec); % assign bins given by pipSecBinVec
            nEventsCell{cc} = accumarray(pipbins,nEvents,[],@(x) {x}); % given cell contains all time bins that have a certain pips/s
            nf_mean{cc} = cellfun(@mean,nEventsCell{cc});
            for bb = 1:length(nEventsCell{cc}) % bootci does not take cell arrays
                nf_ci{cc}(:,bb) = bootci(1000,@mean,nEventsCell{cc}{bb});
            end

            xvec{cc} = pipSecBinCenters(unique(pipbins));
            % stats: is this roi/condition modulated by pip rate?
            p1(rr,cc) = kruskalwallis(nEvents,pipbins,'off'); 
            cc = cc+1;
        end
    end % end D1/D2

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

    sigPipMod = p1(rr,:) < sigThr_pipMod
    sigDistrF1 = squeeze(p2(1,rr,:) < sigThr_audOnly)'
    sigDistrF2 = squeeze(p2(2,rr,:) < sigThr_audOnly)'
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
    pause
end






