
clear
close all
clc


% initialize some other variables:
slugInit


% The observations: 

% Head [m]:
H = [0.72,0.49,0.30,0.20,0.16,0.12];

% measurement times:
TM = [5.0,10.0,20.0,30.0,40.0,50.0];




% simulation times:
TS = [0:0.1:100];

% known parameter values:
D = 60;
Q = 50;


hFig=figure('numbertitle','off','name','manual calibration');
[hEdit1,hEdit2,hButtonGo,hButtonClf,activeSubplot]=disp_ui(hFig);

