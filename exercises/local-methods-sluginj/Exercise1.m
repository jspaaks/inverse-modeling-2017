% Global variables for observed Theta and h values, etc.
%
global H;
global TM;
global Q;
global D;

%
% The unknown/estimated parameters are S=p(1) and T=p(2).
%
%
% The data.  Head is measured to the nearest centimeter.  
%
H=[0.55 0.47 0.30 0.22 0.17 0.14];
TM=[5.0 10.0 20.0 30.0 40.0 50.0];
TS=[0:0.1:60];

%
% Fixed parameter values.
%
D=10;
Q=50;
%
