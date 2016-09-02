% do some p value sanity checks

clear;clc
basedir = '/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone/';
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));
fileID = fopen([basedir 'dirList.txt'],'r');
tmp = textscan(fileID,'%s %d');
dirList = tmp{1};
depthList = tmp{2};

correctN = 1;
s = 0.01;

for ii = 1:length(dirList)
    p = load([basedir dirList{ii} filesep 'pvalues.mat']);
    if correctN
        a1 = s/numel(p.p1);
        a2 = s/numel(p.p2);
    else
        a1 = s/size(p.p1,2);
        a2 = s/(size(p.p2,1)*size(p.p2,3));
    end
    
    subplot(1,2,1);
    p1 = p.p1(p.p1<0.2);
    histogram(p1(:),100)
    axis square tight
    xlim([0 0.2]);
    subplot(1,2,2);
    p2 = p.p2(p.p2<0.2);
    histogram(p2(:),100)
    axis square tight
    xlim([0 0.2]);
    pause
end