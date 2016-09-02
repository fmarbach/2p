function [vec,tr]=dj_trial_extract(dj)
% ao=1;
% sensor=2;
% pulse=3;
% in_zone=4;
% sync=5;
% ch2pulse=6;
% rewards=7;
% trialOn=8;
% licks=9;

vec.sr=500;%hz
vec.stage=dj.params(end,11);
vec.convwin=22;%% ms
vec.convwinwidth=vec.convwin*(vec.sr/1000);%% ms
vec.t=dj.t;

vec.ao=dj.trace(:,1);
vec.sensor=dj.trace(:,2);
vec.pulse=dj.trace(:,3);
vec.in_zone=dj.trace(:,4);
vec.clickDelay=dj.trace(:,5);
vec.reward=dj.trace(:,6);
vec.level=dj.trace(:,7);
vec.sync=dj.trace(:,8);
vec.trialOn=dj.trace(:,9);
vec.ch2_pulse=dj.trace(:,10);
vec.fake_level=dj.trace(:,11);

if size(dj.trace,2)>11
vec.TfType=dj.trace(:,13);
vec.lightType=dj.trace(:,12);
end
vec.fake_level_on=vec.fake_level>0;

vec.rewards_cumsum=cumsum(vec.reward);
vec.block=ones(size(vec.t));%% initiate vec


%%% psiol
try
vec.t_spikes=dj.t_spikes;
end
%%% extract data from param table

%%% extract data from param table
% vec.click_delay=ones(size(vec.ao))*-1;
for ii=1:size(dj.params,1)
    if ii==size(dj.params,1)%% last block
        inx=dj.t>dj.params(ii,1) & dj.t<=max(dj.t);
    else
        inx=dj.t>dj.params(ii,1) & dj.t<=dj.params(ii+1,1);
    end
%     vec.click_delay(inx)=dj.params(ii,12);
    vec.ch2_order(inx,1)=dj.params(ii,14);
    vec.ch2_delay(inx,1)=dj.params(ii,13);
    vec.ch2_block(inx,1)=dj.params(ii,15);
    try%% to aviod error where these field dont exist
    vec.freq1(inx,1)=dj.params(ii,16);
    vec.freq2(inx,1)=dj.params(ii,17);
    end
    vec.block(inx)=ii;
    
       vec.reward_in_block(ii)=sum(vec.reward(inx));%%% check how many reward in the first  block only
    
end

vec.reward_in_block=max(vec.reward_in_block);
vec.uni_gain=unique(dj.params(:,10));
vec.uni_level=unique(dj.params(:,4));
vec.uni_clickdelay=unique(vec.clickDelay);vec.uni_clickdelay(vec.uni_clickdelay==-1)=[];
vec.uni_ch2_delay=unique(dj.params(:,13));
vec.uni_order=unique(dj.params(:,14));

try%% to aviod error where these field dont exist
vec.uni_freq1=unique(vec.freq1);
vec.uni_freq2=unique(vec.freq2);
end


vec.intWin=dj.params(end,6);



vec.uni_clickdelay=unique(vec.clickDelay);vec.uni_clickdelay(vec.uni_clickdelay==-1)=[];
vec.rewards_cum_block=zeros(size(vec.reward));


for ii=1:max(vec.block)
% num_of_rew(ii)=sum(vec.reward(vec.block==ii));
vec.rewards_cum_block(vec.block==ii)=cumsum(vec.reward(vec.block==ii));
end

for ii=1:length(vec.uni_level)
linx=find(dj.params(1:end,4)==vec.uni_level(ii),1,'last');
vec.winwidth(ii,:)=[dj.params(linx,3) dj.params(linx,5)];
end
vec.winLong=dj.params(end,6);

%%% base and the width that were use get the base width (and hypothetical
%%% zeros)
if size(vec.winwidth,1)>1
a=(vec.uni_level(2)-vec.uni_level(1)).\(vec.winwidth(2,2)-vec.winwidth(1,2));
vec.winWide=vec.winwidth(1,2)-vec.uni_level(1)*a;
else
 vec.winWide=vec.winwidth(1,2);
end


%%% get reward in block times
vec.rwd.inx=find(vec.reward);
for rr=1:length(vec.rwd.inx)
    vec.rwd.block(rr,1)=vec.block(vec.rwd.inx(rr));
    vec.rwd.time_in_block(rr,1)=vec.t(vec.rwd.inx(rr))- dj.params(vec.rwd.block(rr),1);
    vec.rwd.cum_num(rr,1)=vec.rewards_cum_block(vec.rwd.inx(rr));
    if vec.rwd.cum_num(rr,1)==1
        vec.rwd.iri(rr,1)=vec.t(vec.rwd.inx(rr))-dj.params(vec.rwd.block(rr),1);
        
    else
        vec.rwd.iri(rr,1)=vec.t(vec.rwd.inx(rr))-vec.t(vec.rwd.inx(rr-1));
    end
end


%%% get block data

vec.blck.level=dj.params(1:end-1,4);
vec.blck.dur=diff(dj.params(:,1));
vec.blck.delay=dj.params(1:end-1,12);
vec.blck.num=[1:size(vec.blck.level,1)]';
vec.blck.ch2_order=dj.params(1:end-1,14);
for ii=1:size(vec.blck.level,1)
    tinx=find(vec.t>dj.params(ii,1),1):find(vec.t<dj.params(ii+1,1),1,'last');
   vec.blck.reward(ii,1)=sum(vec.reward(tinx));
end


%%% get transfer functions
offtarget=linspace(6,24,2);
ontarget=linspace(48,96,2);
vec.tf=repmat([offtarget ontarget fliplr(ontarget) fliplr(offtarget)],length(vec.uni_level),1)';
for ii=1:length(vec.uni_level)
    vec.tf_axis(:,ii)=[0,vec.uni_level(ii)+vec.winwidth(ii,1), ...
        vec.uni_level(ii)+vec.winwidth(ii,1),vec.uni_level(ii), ...
        vec.uni_level(ii),vec.uni_level(ii)+vec.winwidth(ii,2) ...
        vec.uni_level(ii)+vec.winwidth(ii,2),5];
end

    %%% rate states
vec.rate_state=zeros(size(vec.ao));
pulse_locs=find(vec.pulse);
pulse_rate=diff([NaN;vec.t(pulse_locs)]);
pulse_rate_state=pulse_rate;
pulse_rate_state(pulse_rate_state>.1)=0;
for ii=2:length(pulse_locs)
    vec.rate_state(pulse_locs(ii-1)+1:pulse_locs(ii))=pulse_rate_state(ii);
    vec.rate(pulse_locs(ii-1)+1:pulse_locs(ii))=pulse_rate(ii);
end
vec.rate_state(vec.rate_state>.041)=1;
vec.rate_state(vec.rate_state<=.041 & vec.rate_state>0)=2;
vec.rate(vec.rate>1/6)=1/3;
%%%% ch2 %%%%%%
vec.ch2_rate_state=zeros(size(vec.ao));%%%% ch1 
vec.ch2_rate=zeros(size(vec.ao));%%%% ch1 
pulse_locs=find(vec.ch2_pulse);
pulse_rate=diff([NaN;vec.t(pulse_locs)]);
pulse_rate_state=pulse_rate;
pulse_rate_state(pulse_rate_state>.1)=0;
for ii=2:length(pulse_locs)
    vec.ch2_rate_state(pulse_locs(ii-1)+1:pulse_locs(ii))=pulse_rate_state(ii);
    vec.ch2_rate(pulse_locs(ii-1)+1:pulse_locs(ii))=pulse_rate(ii);
end
vec.ch2_rate_state(vec.ch2_rate_state>.041)=1;
vec.ch2_rate_state(vec.ch2_rate_state<=.041 & vec.ch2_rate_state>0)=2;
vec.ch2_rate(vec.ch2_rate>1/6)=1/3;

%%% reward rate
% sessionInMinutes=floor(vec.t(end)/60);
% vec.rewardPerMinuteRate=sum(reshape(vec.reward(1:vec.sr*sessionInMinutes*60),vec.sr*60,sessionInMinutes));
% trsh=5;%reward per minutes;
% % vec.miceIsWorking=vec.t<(find(vec.rewardPerMinuteRate>trsh,1,'last')*60);

%% make the freq vec
vec.f1=zeros(size(vec.pulse));
vec.f2=zeros(size(vec.pulse));
vec.f1(vec.ch2_order==12)=vec.pulse(vec.ch2_order==12);
vec.f1(vec.ch2_order==21)=vec.ch2_pulse(vec.ch2_order==21);
vec.f2(vec.ch2_order==21)=vec.pulse(vec.ch2_order==21);
vec.f2(vec.ch2_order==12)=vec.ch2_pulse(vec.ch2_order==12);
%% trials
trialUp=find(diff(vec.trialOn)==1)+1;
trialDw=find(diff(vec.trialOn)==-1)+1;
if trialDw(1)<trialUp(1)%% make sure to have only full trilas
    trialDw(1)=[];
end
if trialDw(end)<trialUp(end)
    trialUp(end)=[];
end

if trialUp(1)<1500%% delete first trail if on in the first 3 seconds
    trialUp(1)=[];
    trialDw(1)=[];
end

for ii=1:length(trialUp)
   tr.t{ii,1}=vec.t(trialUp(ii):trialDw(ii))-vec.t(trialUp(ii));
   tr.ao{ii,1}=vec.ao(trialUp(ii):trialDw(ii)); 
   
   trueCross=find(vec.ao(trialUp(ii)-199:trialUp(ii))>4.7,1,'last');
   if ~any(trueCross)
       trueCross=200;
   end
   corrOnset=trialUp(ii)-(200-trueCross);
   tr.aoCr{ii,1}=vec.ao(corrOnset:trialDw(ii)); 
   tr.tCr{ii,1}=vec.t(corrOnset:trialDw(ii))-vec.t(trialUp(ii));
if size(dj.trace,2)>11
   tr.TfType(ii,1)=vec.TfType(trialUp(ii)); 
   tr.lightType(ii,1)=vec.lightType(trialUp(ii)); 
end
   tr.pulse{ii,1}=vec.pulse(trialUp(ii):trialDw(ii)); 
   tr.ch2_pulse{ii,1}=vec.ch2_pulse(trialUp(ii):trialDw(ii)); 
   tr.block(ii,1)=vec.block(trialUp(ii));
   tr.rewards_cum_block(ii,1)=vec.rewards_cum_block(trialUp(ii));
   tr.outcome(ii,1)=sum(vec.reward([-10:10]+trialDw(ii)));
   tr.level(ii,1)=vec.level(trialUp(ii));
   tr.fake_level(ii,1)=vec.fake_level(trialUp(ii));
   tr.clickDelay(ii,1)=vec.clickDelay(trialUp(ii));
   tr.inZone{ii,1}=vec.in_zone(trialUp(ii):trialDw(ii)); 
   tr.t0(ii,1)=vec.t(trialUp(ii));
   tr.aoWarp(ii,:)=interp1(tr.tCr{ii},tr.aoCr{ii},linspace(tr.tCr{ii}(1),tr.t{ii}(end),100));
   tr.ch2Order(ii,1)=vec.ch2_order(trialUp(ii));
   if isfield(dj,'spikes')
   tr.spikes{ii,1}=dj.spikes(trialUp(ii)*2:trialDw(ii)*2,:);
   end
end
tr.fake_level_on=tr.fake_level>0;
tr.num=[1:length(tr.outcome)]';
tr.rewards_cum_block=tr.rewards_cum_block+1;
tr.rt=cellfun(@length,tr.t)*1/vec.sr;


% vec.miceIsWorking=zeros(size(vec.t));
% if any(find(vec.rewardPerMinuteRate>trsh,1,'last')*60)
% vec.miceIsWorking(vec.t<tr.t0(find(tr.outcome==1,1,'last')))=1;%% from beggingin to last succesful trial
% else
%   vec.miceIsWorking=zeros(size(vec.t))  ;
% end
