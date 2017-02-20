%
% Reset to a known state.
%
clear
rand('state',0);
randn('state',0);
%
% Global variables for observed Theta and h values, etc.
%
global H
global TM
global TS
global Q
global D

%
% The unknown/estimated parameters are S=p(1) and T=p(2).
%
%
