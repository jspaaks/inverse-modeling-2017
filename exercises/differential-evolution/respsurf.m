% differential evolution script

% restore the default MATLAB path in case it was changed:
restoredefaultpath

% add the folder to the MATLAB path that contains the 
% data needed for some of the benchmark functions:
addpath('.\benchmark-data')

% clear any old variables:
clear

% close any old figures:
close all

% clen up the command window, so we don't get 
% confused by any old error messages:
clc


% the benchmark function needs a few variables:
global initial_flag; 
initial_flag = 0;           

% this next command sets the flag, indicating which function is used:
funcFlag = 6;

parVec = [0.2,4.7]

objScore = benchmark_func(parVec,funcFlag)


x = [65:0.1:80];
y = [35:0.1:45];

nx = numel(x);
ny = numel(y);

respSurf = repmat(NaN,[ny,nx]);
for iy = 1:ny
    for ix=1:nx
        parVec = [x(ix),y(iy)];
        respSurf(iy,ix) = benchmark_func(parVec,funcFlag);
    end
end

figure
imagesc(x,y,respSurf)
axis image
set(gca,'ydir','normal')
title('objScore')
xlabel('x')
ylabel('y')
colorbar

set(gcf,'paperpositionmode','auto','inverthardcopy','off')
print('respsurf.eps','-depsc2','-r300','-loose')
print('respsurf.png','-dpng','-r300')
