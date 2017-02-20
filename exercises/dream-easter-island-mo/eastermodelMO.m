function OUT = eastermodelMO(parVec,Extra)
% easter island 
% resources and population model in closed system


parPopGrowth = parVec(1);        % population growth rate = birth rate
parDeathRate = parVec(2);    % population death rate
parCarCap = parVec(3);       % carrying capacity
parConsResponse = parVec(4); % consumption response
parRelResGrowth = parVec(5);        % relative resource growth
parFit = parVec(6);       % fitting parameter

tMeas = Extra.tMeas;
PopMeas = Extra.PopMeas;
ResMeas = Extra.ResMeas;


tStart = Extra.tStart;
tEnd = Extra.tEnd;
dt = Extra.dt;

Pop = 20;
Resource = 160;

S = repmat(NaN,[numel(tStart:dt:tEnd)-1,3]);
k = 0;


t = tStart;

while t<tEnd
    
    consRate = parFit*Pop*Resource;
    consRatePerCapita = consRate/Pop;
    consProfit = consRatePerCapita*parConsResponse;
    resourceChangeRate = parRelResGrowth*Resource*(1 - Resource/parCarCap);
    Resource = Resource + resourceChangeRate*dt - consRate*dt;
    
    popGrowth = (parPopGrowth + consProfit)*Pop;
    popDeaths = (parDeathRate - consProfit)*Pop;
    Pop = Pop + popGrowth*dt - popDeaths*dt;
    
    t = t + dt; 
    k = k + 1;
    S(k,1:3) = [t,Pop,Resource];
    if isnan(Pop)
        break
    end
end

if isnan(S(end,2))
    clear S
    S(:,1) = tMeas;
    S(:,2:3) = -999;
end

tSim = S(:,1);
PopSim = S(:,2);
ResSim = S(:,3);


% Interpolate the simulated Population and Resources to the times
% at which measurements are available against which to compare the
% simulations:
PopSimI = interp1(tSim,PopSim,tMeas);
ResSimI = interp1(tSim,ResSim,tMeas);


OUT = PopSimI;


if ~Extra.dontPlot

    % every ten seconds, there is a time window of 1 second, during which the
    % simulation results are visualized. The other 9 seconds, nothing is
    % visualized to speed up the optimization.
    c = clock;
    plotInterval = 10; %s
    if mod(floor(c(6)),plotInterval)==0

        measStyle = {'marker','o','markersize',3,'linestyle','none'};
        figure(123)
        clf
        subplot(2,1,1)
        hResObs = plot(tMeas,ResMeas,'markeredgecolor',[0,0.5,0],...
            'markerfacecolor',[0,0.5,0],measStyle{:});
        hold on
        hResSim = plot(tSim,ResSim,'linestyle','-','color',[0,0.5,0]);
        set(gca,'xlim',[0,2000],'ylim',[0,250])

        xlabel('Time [years]')
        ylabel('Biomass')
        subplot(2,1,2)
        hPopObs = plot(tMeas,PopMeas,'markeredgecolor',[1,0,1],...
            'markerfacecolor',[1,0,1],measStyle{:});
        hold on
        hPopSim = plot(tSim,PopSim,'linestyle','-','color',[1,0,1]);
        set(gca,'xlim',[0,2000],'ylim',[0,15000])
        
        xlabel('Time [years]')
        ylabel('Population')


        drawnow
    end
end
