function activeSubplot=plot_mancal_callback(hEdit1,hEdit2,hEdit3,hEdit4,activeSubplot)



s1 = get(hEdit1,'string');
s2 = get(hEdit2,'string');
s3 = get(hEdit3,'string');
s4 = get(hEdit4,'string');


parVec(1)=str2double(s1);
parVec(2)=str2double(s2);
parVec(3)=str2double(s3);
parVec(4)=str2double(s4);

if isnan(parVec(1))
    blink(hEdit1)
    return
end

if isnan(parVec(2))
    blink(hEdit2)
    return
end

if isnan(parVec(3))
    blink(hEdit3)
    return
end

if isnan(parVec(4))
    blink(hEdit4)
    return
end




% if activeSubplot==4
%     activeSubplot=activeSubplot+1;
% end



% load('bound_211.mat')
% load('stor_211.mat')

% TimeTab = bound_211(:,1);
% PrecTab = bound_211(:,2);
% PEvapTab = bound_211(:,3);
% 
% storageMeasured = stor_211(:,2);

load('bound_ideal.mat')
load('stor_ideal.mat')

TimeTab = bound_ideal(:,1);
PrecTab = bound_ideal(:,2);
PEvapTab = bound_ideal(:,3);

storageMeasured = stor_ideal(:,2);

if ~isequal(stor_ideal(:,1),bound_ideal(:,1))
    error('input data inconsistent.')
end


precision = 1e-6;
TMP=unique(round(((TimeTab(2:end)-TimeTab(1:end-1))/precision))*precision);
if numel(TMP)~=1
    error('variable time step in input data.')
else
    dt = TMP;
end


StartTime = TimeTab(1);
EndTime = TimeTab(end);
Time = StartTime;
Storage(1) = 0;%0.128;

clear bound_211

k=1;

a = parVec(1);
b = parVec(2);
c = parVec(3);
d = parVec(4);

 
while Time<EndTime

  PEvapR = PEvapTab(k);
  PrecR = PrecTab(k);
  St = (Storage(k)/dt+a*PrecR+b*c)/(1/dt+b+d*PEvapR/c);

  if (St<c)
      St = (Storage(k)/dt+a*PrecR)/(1/dt+d*PEvapR/c);
      if St>c
          St=c;
      end
  end
    
   k = k + 1;
   Storage(k) = St;
   Time = Time + dt;
   Storsave(k-1) = (Storage(k)-Storage(k-1))./2+Storage(k-1);
        
end
% Storsave(k) = Storsave(k-1);
% S = Storsave';

% technically it's not allowed to export this variable rather than
% Storsave (due to implicit integration):
S = Storage;


subplot(2,3,activeSubplot)

yMax = 3.5;

h=plot(StartTime:dt:EndTime,yMax-0.9*scale(PrecTab)*yMax,'-c',...
     StartTime:dt:EndTime,0.9*scale(PEvapTab)*yMax,'-r',...
     StartTime:dt:EndTime,storageMeasured,'mo',...
     StartTime:dt:EndTime,Storage,'-b.');
 
set(h(3),'markersize',3,'markerfacecolor','w') 
set(gca,'ylim',[0,yMax],'xlim',[StartTime,EndTime])
xlabel('Time [days]')
ylabel('Canopy storage');
title(['a = ',num2str(parVec(1)),' ; b = ',num2str(parVec(2)),...
    ' ; c = ',num2str(parVec(3)),' ; d = ',num2str(parVec(4))])


if activeSubplot==1
%     subplot(2,3,4)
% 
%     yMax = 3.5;
% 
%     h=plot(StartTime:dt:EndTime,yMax-0.9*scale(PrecTab)*yMax,'-c',...
%          StartTime:dt:EndTime,0.9*scale(PEvapTab)*yMax,'-r',...
%          StartTime:dt:EndTime,storageMeasured,'mo',...
%          StartTime:dt:EndTime,S,'-b.');
% 
%     set(h(3),'markersize',3,'markerfacecolor','w')
%     set(gca,'xlim',[StartTime,EndTime]-StartTime-100)
%     activeSubplot= activeSubplot+1;
    hLeg=legend(h,'scaled precipitation','scaled evaporation',...
        'measured storage','simulated storage','location','NorthWest');
    set(hLeg,'color','none','edgecolor','w','fontsize',8)
%     axis off
% 
end





if activeSubplot<8
    activeSubplot=activeSubplot+1;
end



function blink(hEdit)

for u=1:4
    set(hEdit, 'BackgroundColor','r')
    drawnow
    pause(0.25)
    set(hEdit, 'BackgroundColor','w')
    pause(0.25)
    drawnow
end

function Z = scale(Y)


y0 = min(Y(:));
y1 = max(Y(:));


Z = (Y-y0)./(y1-y0);
