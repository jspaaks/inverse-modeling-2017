function output = interceptionmodel(parVec,Extra)
% Single layer interception model
% rvw: set met randvoorwaarden
% Column 1: Time [day]
% Column 2: Cumulative precipitation [mm]
% Column 3: Cumulative potential evaporation [mm]

a = parVec(1);
b = parVec(2);
c = parVec(3);
d = parVec(4);


TimeTab = Extra.bound(:,1);
precision = 1e-6;
TMP=unique(round(((TimeTab(2:end)-TimeTab(1:end-1))/precision))*precision);
if numel(TMP)~=1
    error('variable time step in input data.')
else
    dt = TMP;
end


PrecTab = Extra.bound(:,2);
PEvapTab = Extra.bound(:,3);

Storage(1) = Extra.stor(1,2);


k=1;
dt=TimeTab(2:end)-TimeTab(1:end-1);
StartTime = TimeTab(1);
EndTime = TimeTab(end); 

Time = StartTime;
%for t=StartTime:dt:(EndTime)
while Time< EndTime %%%%%%%%
  PEvapR = PEvapTab(k);
  PrecR = PrecTab(k);
  St = (Storage(k)/dt(k)+a*PrecR+b*c)/(1/dt(k)+b+d*PEvapR/c);

  if (St<c)
      St = (Storage(k)/dt(k)+a*PrecR)/(1/dt(k)+d*PEvapR/c);
      if St>c
          St=c;
      end
  end
    
   k=k+1;
   Storage(k) = St;
   Time = Time + dt;
   Storsave(k-1) = (Storage(k)-Storage(k-1))./2+Storage(k-1);
        
end
% Storsave(k) = Storsave(k-1);
% S = Storsave';

% technically it's not allowed to export this variable rather than
% Storsave (due to implicit integration):
S = Storage;

output = Storage';




