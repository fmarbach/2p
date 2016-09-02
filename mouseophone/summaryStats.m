% load pvalues.mat from different sessions and combine results

clear;close all;clc
basedir = '/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone/';
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));

% read text file with experiment info
fileID = fopen([basedir 'dirList.txt'],'r');
tmp = textscan(fileID,'%s %d');
dirList = tmp{1};
depthList = tmp{2};

% read mat file with pip rate stats
load([basedir 'pipDiff.mat']);

correctN = 1;
s = 0.01;

for ii = 1:length(dirList)
    p = load([basedir dirList{ii} filesep 'pvalues.mat']);
    t = load([basedir dirList{ii} filesep 'trial.mat']);
    if correctN
        a1 = s/numel(p.p1);
        a2 = s/numel(p.p2);
    else
        a1 = s/size(p.p1,2);
        a2 = s/(size(p.p2,1)*size(p.p2,3));
    end
    
    % p1: modulated by pip rate?
    % p2: distractor different from feedback? (for given pip rate)
    p1H = p.p1 < a1;
    p1H_F1 = sum(p1H(:,1:2),2); % freq1 distractor and feedback
    p1H_F2 = sum(p1H(:,3:4),2); % freq2 distractor and feedback
    
    p2H = p.p2 < a2;
    p2H_F1 = sum(squeeze(p2H(1,:,:)),2); % freq1
    p2H_F2 = sum(squeeze(p2H(2,:,:)),2); % freq2
    
    nCells = length(p2H_F1);
    aud = p1H_F1 | p1H_F2; % aud
    feed = (p1H_F1 & p2H_F1) | (p1H_F2 & p2H_F2); % (aud & feed) for either freq
    nomod = ~aud; 
    audC = sum(aud);
    feedC = sum(feed);
    nomodC = sum(nomod);
    allpie(ii,:) = 100*[nomodC audC-feedC feedC]/nCells;
    h(ii) = subplot(3,ceil(length(dirList)/3),ii);
    pie(round(allpie(ii,:)),{['nomod ' num2str(round(allpie(ii,1))) '%'],...
                             ['aud ' num2str(round(allpie(ii,2))) '%'],...
                             ['feed ' num2str(round(allpie(ii,3))) '%']});
    
    pS(ii) = t.pSuccess;
    
end

figure;
allpieMean = round(mean(allpie,1));
pie(allpieMean,{['nomod ' num2str(allpieMean(1)) '%'],...
                    ['aud ' num2str(allpieMean(2)) '%'],...
                    ['feed ' num2str(allpieMean(3)) '%']});
    
figure;
lab = {'nomod','aud','feed'};
for ii = 1:3
    
    subplot(2,3,ii);
    line([0 750],[0 100]);
    plot(depthList,allpie(:,ii),'k*'); hold on
    axis tight square;
    title(['depth vs. ' lab{ii}]);
    
    subplot(2,3,ii+3);  
    plot(pS,allpie(:,ii),'k*'); hold on
    axis tight square;
    title(['% success vs. ' lab{ii}]);
end
    

figure;
plot(pipDiff,allpie(:,3),'k*')
formatFigure(gcf, gca, 12, 2, 2, 'Pip rate discrepancy [pips/s]', '% feedback modulated cells', ...
             0, 0, 1, 0, 0, 1, 1)    
axis square tight