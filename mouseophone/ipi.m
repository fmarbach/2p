% plot some pip statistics feedback vs. distractor

clear;clc
basedir=uigetdir('/Users/fmarbach/Documents/WSBS/zadorlab/DATA/IMAGING/mouseophone');
addpath(genpath('/Users/fmarbach/Dropbox/CODE/2p_Margaret/mouseophone/'));
[~,sessionName]=fileparts(basedir);
load([basedir filesep sessionName '.mat'])
[vec,tr] = dj_trial_extract(dj);

%%
close all;
thr = 0.5;
w = 5;

% ensure that we start trial off and end trial off
tmp = vec.trialOn(find(vec.trialOn==0,1,'first'):find(vec.trialOn==0,1,'last'));
trialOff = find(diff(tmp)==-1);
trialOn = find(diff(tmp)==1);

r1 = [];
r2 = [];
rF = [];
rD = [];
ipi1 = {};
ipi2 = {};
block = [];
ind = 0;
for tt = 1:length(trialOn)
    s = trialOn(tt):trialOff(tt);
    if length(s)/500 > thr
        ind = ind+1;
        r1(ind) = sum(vec.f1(s))/(length(s)/500);
        r2(ind) = sum(vec.f2(s))/(length(s)/500);
        rF(ind) = sum(vec.pulse(s))/(length(s)/500);
        rD(ind) = sum(vec.ch2_pulse(s))/(length(s)/500);
%         ipi1{ind} = diff(find(vec.f1(s)==1));
%         ipi2{ind} = diff(find(vec.f2(s)==1));
        block(ind) = mode(vec.ch2_order(s)==12);
    end
end

f2 = figure(2);
%plot(r1,r2,'ko'); hold on
%line([0 max([r1 r2])],[0 max([r1 r2])]);
plot(r1,'ko'); hold on;
plot(r2,'b*');
plot(medfilt1(r1,w),'k-','LineWidth',2);
plot(medfilt1(r2,w),'b-','LineWidth',2);
plot(block*max([r1 r2]),'r-');
formatFigure(gcf, gca, 12, 2, 2, 'trial number', 'pips/sec', ...
             0, 0, 1, 0, 0, 1, 1)
axis tight;

f3 = figure(3);
xvec = (1:length(vec.f1))/500;
mr1 = smooth(vec.f1,5000);
mr2 = smooth(vec.f2,5000);
plot(xvec,mr1); hold on
plot(xvec,mr2);
plot(xvec,max([mr1; mr2])*(vec.ch2_order==12))
    
f4 = figure(4);
%plot(r1,r2,'ko'); hold on
%line([0 max([r1 r2])],[0 max([r1 r2])]);
plot(r1-r2,'ko'); hold on;
plot(medfilt1(r1-r2,w),'k-','LineWidth',2);
plot(block*max([r1-r2]),'r-');
formatFigure(gcf, gca, 12, 2, 2, 'trial number', 'pips/sec', ...
             0, 0, 1, 0, 0, 1, 1)
axis tight;

f5 = figure(5);
histogram(rF-rD,100)
mean(rF-rD)
    