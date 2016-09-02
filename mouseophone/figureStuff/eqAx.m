function eqAx(correct,handle_list)
if nargin < 1
    correct='both';
end
H=get(gcf,'Children');
if nargin==2
    H=handle_list;
end

hValid=zeros(size(H));
for  ii=1:length(H)
    if strcmp(H(ii).Type,'legend')
        hValid(ii)=1;
    end
end
H(hValid==1)=[];


for h=1:length(H)
    val(h)=sum(get(H(h),'YLim')=='on')~=2;
    ylim(h,:)=get(H(h),'YLim');
    xlim(h,:)=get(H(h),'XLim');
        zlim(h,:)=get(H(h),'ZLim');
end
ymin=min(ylim(val,1));ymax=max(ylim(val,2));
xmin=min(xlim(val,1));xmax=max(xlim(val,2));
zmin=min(zlim(val,1));zmax=max(zlim(val,2));
for h=find(val)
    if correct=='both' | correct=='y';set(H(h),'YLim',[ymin ymax]);end
    if correct=='both' | correct=='x';set(H(h),'XLim',[xmin xmax]);end
    if correct=='both' | correct=='z';set(H(h),'ZLim',[zmin zmax]);end

end


