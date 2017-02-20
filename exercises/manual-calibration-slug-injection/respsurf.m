
clear
clear global
close all
clc


% initialize some other variables:
slugInit


% The observations: 

% measurement times:
TM = [5.0,10.0,20.0,30.0,40.0,50.0];
% measured head [m]:
HM = [0.72,0.49,0.30,0.20,0.16,0.12];


% simulation times:
TS = [0.1:0.1:100];

% known parameter values:
D = 60;
Q = 50;



S = linspace(0.0001,0.01,151);
T = linspace(0.01,2.0,101);

nS = numel(S);
nT = numel(T);

% pre-allocate respSurf for speed:
respSurf = repmat(NaN,[nT,nS]);

for r=1:nT
    for c=1:nS    
        
        parVec = [S(c),T(r)];
        
        % simulated head:
        HS = sluginj(parVec,TS);
        
        % interpolated HS to get the simulation at the same 
        % time that the measurement was taken:
        HSI = interp1(TS,HS,TM);
        
        % objective function is sum of squared residuals:
        respSurf(r,c) = sum((HSI-HM).^2);
        
        
        % figure(123)
        % plot(TS,HS,'-k',TM,HSI,'ok',TM,HM,'.m')
        % set(gca,'ylim',[0,1])
        % drawnow
        
        
    end
end



figure
imagesc(S,T,respSurf,[0,0.5])
hold on
plot(0.00207+[-1,-1]*0.00012,[0,2],'--',...
     0.00207+[1,1]*0.00012,[0,2],'--','linewidth',1,'color',0.9*[1,1,1]);
plot([0,0.1],0.585+[1,1]*0.029,'--',...
     [0,0.1],0.585+[-1,-1]*0.029,'--','linewidth',1,'color',0.9*[1,1,1]);
set(gca,'ydir','normal')
colorbar
title('response surface for slug injection (SSR)')
xlabel('S')
ylabel('T')


set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('respsurf-sluginj.eps','-depsc2','-r300','-loose')
print('respsurf-sluginj.png','-dpng','-r300')


