clear 
close all
clc


for flag = 1:5

    switch flag
        case 1
            noiseFactor = 0.01;
            matFileName = 'pop-res-data-ref.mat';
            t700change = false;
            highDeathRate = false;
        case 2
            noiseFactor = 0.20;
            matFileName = 'pop-res-data-ref-noise.mat';
            t700change = false;
            highDeathRate = false;
        case 3
            noiseFactor = 0.01;
            matFileName = 'pop-res-data-variable-deathrate.mat';
            t700change = true;
            highDeathRate = false;
        case 4
            noiseFactor = 0.20;
            matFileName = 'pop-res-data-variable-deathrate-noise.mat';
            t700change = true;
            highDeathRate = false;            
        case 5
            noiseFactor = 0.01;
            matFileName = 'pop-res-data-high-deathrate.mat';
            t700change = false;
            highDeathRate = true;            
    end

    % easter island 
    % resources and population model in closed system


    parPopGrowth = 0.02;        % population growth rate = birth rate
    if highDeathRate
        parDeathRate = 0.035;   % population death rate
    else
        parDeathRate = 0.03;    % population death rate        
    end
    parCarCap = 160;            % carrying capacity
    parConsResponse = 150;      % consumption response
    parRelResGrowth = 0.004;    % relative resource growth
    parFit = 4.2e-7;            % fitting parameter



    t = 0;
    tStart = 0;
    tEnd = 2000;
    dt = 1.0;

    Pop = 20;
    Resource = 160;

    S = repmat(NaN,[numel(tStart:dt:tEnd)-1,3]);
    k = 0;

    while t<tEnd
        
        
        if t700change && t>700
            parDeathRate = 0.035;
        end


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
    ResourceSim = S(:,3);


    randn('state',0)
    PopMeasNoisy = PopSim.*(1+noiseFactor*randn(size(PopSim)));
    ResMeasNoisy = ResourceSim.*(1+noiseFactor*randn(size(ResourceSim)));

    
    rand('twister',0)
    tmp = randperm(2000);
    sel = sort(tmp(1:40));

    
    tMeas = S(sel,1);
    PopMeas = PopMeasNoisy(sel,1);
    ResourceMeas = ResMeasNoisy(sel,1);


    measStyle = {'marker','o','markersize',3,'linestyle','none'};
    figure(123)
    clf
    subplot(2,1,1)
    hResSim = plot(tSim,ResourceSim,'color',[0,0.5,0],...
        'linestyle','-');
    hold on
    hResMeas = plot(tMeas,ResourceMeas,'markeredgecolor',[0,0.5,0],...
        'markerfacecolor',[0,0.5,0],measStyle{:},'linestyle','none');

    xlabel('Time [years]')
    ylabel('Biomass')

    subplot(2,1,2)
    hPopSim = plot(tSim,PopSim,'color',[1,0,1],...
        'linestyle','-');
    hold on
    hPopMeas = plot(tMeas,PopMeas,'markeredgecolor',[1,0,1],...
        'markerfacecolor',[1,0,1],measStyle{:},'linestyle','none');
    xlabel('Time [years]')
    ylabel('Population')
    drawnow


    save(matFileName,'tMeas','PopMeas','ResourceMeas')

end



