clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

[~,sessionName]=fileparts(basedir);
load([basedir filesep sessionName '.mat'])
[vec,tr] = dj_trial_extract(dj);

longThreshold = 1.1;

nTrials = length(tr.lightType);
light = tr.lightType == 1;
long = tr.rt > longThreshold;

p = sum(tr.outcome) / nTrials;
pNoLight = sum(tr.outcome(~light)) / sum(~light);
pLight = sum(tr.outcome(light)) / sum(light);

p_L = sum(tr.outcome(long)) / sum(long);
pNoLight_L = sum(tr.outcome(~light & long)) / sum(~light & long);
pLight_L = sum(tr.outcome(light & long)) / sum(light & long);

p_S = sum(tr.outcome(~long)) / sum(~long);
pNoLight_S = sum(tr.outcome(~light & ~long)) / sum(~light & ~long);
pLight_S = sum(tr.outcome(light & ~long)) / sum(light & ~long);

% durLight = tr.rt(light);
% durNoLight = tr.rt(~light);
% l = histcounts(durLight,100) / length(durLight);
% nl = histcounts(durNoLight,100) / length(durNoLight);
% plot(l); hold on
% plot(nl);

% plot bars
h = bar([p p_L p_S; pLight pLight_L pLight_S; pNoLight pNoLight_L pNoLight_S]);
legend('all','long','short')

%% plot trajectories
thr = 1;
levels = unique(tr.level);
for jj = 1:length(levels)
    t{jj,1} = find(tr.rt>thr & tr.level==jj & tr.lightType & tr.outcome==1);
    t{jj,2} = find(tr.rt>thr & tr.level==jj & ~tr.lightType & tr.outcome==1);
end

for jj = 1:length(levels)
    for ii = 1:2 % light vs no light
        figure(1);
        h = subplot(length(levels),2,2*jj-2+ii); hold on
        len = cellfun(@length,tr.ao(t{jj,ii}));
        pad = NaN(1,max(len));
        traceAvg = zeros(max(len),length(len));
        for tt = 1:length(len)
            trace = tr.ao{t{jj,ii}(tt)};
            tmp = pad;
            tmp(end-length(trace)+1:end) = trace;
            traceAvg(:,tt) = tmp;
            
            %traceFlip = flipud(trace);
            %plot(traceFlip)
            plot(tmp,'Color',0.9*[1 1 1])
        end
        plot(nanmean(traceAvg,2),'k-')
        axis tight
        
        figure(2);
        subplot(length(levels),1,jj);
        plot(nanmean(traceAvg,2)); hold on;
    end
end








