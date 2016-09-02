function [clims]=eqAxC(H,eqlowhigh)
if nargin==0 
H=get(gcf,'Children');
end
if nargin<=1
    eqlowhigh=0;
end
%%% check is colorbar
hValid=zeros(size(H));
for h=1:length(H)
    if strcmp(H(h).Type,'colorbar')
        hValid(h)=1;
    end
end
H(hValid==1)=[];

for h=1:length(H)
    val(h)=sum(get(H(h),'CLim')=='on')~=2;
    clim(h,:)=get(H(h),'CLim');
end
cmin=min(clim(val,1));cmax=max(clim(val,2));
if eqlowhigh
    newlim=max(abs([cmin cmax]));
    cmin=-newlim;
    cmax=newlim;
end

for h=find(val)
    set(H(h),'CLim',[cmin cmax]);
end
clims=[cmin cmax];

