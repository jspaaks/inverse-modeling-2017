function output = hymod(parVec,Extra) 
% Runs the HYMOD model

cmax = parVec(1);
bexp = parVec(2);
fQuickFlow = parVec(3);
Rs = parVec(4);
Rq = parVec(5);

dailyDischarge = Extra.dailyDischarge;
dailyPotEvapTrans = Extra.dailyPotEvapTrans;
dailyPrecip = Extra.dailyPrecip;
iStart = Extra.iStart;
iEnd = Extra.iEnd;
numTime = Extra.numTime;

% HYMOD PROGRAM IS SIMPLE RAINFALL RUNOFF MODEL
% START PROGRAMMING LOOP WITH DETERMINING RAINFALL - RUNOFF AMOUNTS

convFactor = 22.5;

wu = 65;

t = iStart-wu;
tEnd = iEnd;

x_loss = 0.0;

% Initialize slow tank state
x_slow = 2.3503/(Rs*convFactor);

% Initialize state(s) of quick tank(s)
x_quick(1:3,1) = 0;

outflow = [];


while t < tEnd+1
    
   Pval = dailyPrecip(t,1);
   PETval = dailyPotEvapTrans(t,1);
   
   % Compute excess precipitation and evaporation
   [UT1,UT2,x_loss] = excess(x_loss,cmax,bexp,Pval,PETval);
   
   % Partition UT1 and UT2 into quick and slow flow component
   UQ = fQuickFlow*UT2 + UT1; 
   US = (1-fQuickFlow)*UT2;
   
   % Route slow flow component with single linear reservoir
   inflow = US; 
   [x_slow,outflow] = linres(x_slow,inflow,outflow,Rs); 
   QS = outflow;
   
   % Route quick flow component with linear reservoirs
   inflow = UQ; 
   k = 1; 
   while k < 4
      [x_quick(k),outflow] = linres(x_quick(k),inflow,outflow,Rq); 
      inflow = outflow; 
      k = k+1;
   end
   
   % Compute total flow for timestep
   tempOutput(t,1) = (QS + outflow)*convFactor;
   t = t+1;   
end

output=tempOutput(iStart:iEnd);


SimRR = [numTime(iStart:iEnd),output];







% every ten seconds, there is a time window of 1 second, during which the
% simulation results are visualized. The other 9 seconds, nothing is
% visualized to speed up the optimization.
c = clock;
plotInterval = 10; %s
if mod(floor(c(6)),plotInterval)==0

    figure(123)
    subplot(2,1,1)
    plot(numTime(iStart:iEnd),dailyPotEvapTrans(iStart:iEnd),'-r',...
         numTime(iStart:iEnd),0.1*dailyPrecip(iStart:iEnd),'-c' )
    hLeg = legend('dailyPotEvapTrans [mm\cdot{}day^{-1}]','dailyPrecip [10\cdot{}mm\cdotday^{-1}]');
    set(hLeg,'fontsize',8,'color','none','xcolor','w','ycolor','w','location','northwest')
    
    title('Leaf River watershed observations')    
    datetick('x')

    subplot(2,1,2)
    h=plot(numTime(iStart:iEnd),dailyDischarge(iStart:iEnd),'om',...
           SimRR(:,1),SimRR(:,2),'-b');
    set(gca,'ylim',[0,420])   
    set(h(1),'markersize',3,'markerfacecolor','m')
    ylabel('dailyDischarge [m^3\cdot{}s^{-1}]')
    datetick('x')
    drawnow
end