%
% fvec=fun(p)
%
% Computes the differences between forward model prediction and data,
% normalized by the standard deviation for the slug test example.
%
function fvec=fun(p)
%
% global variables.
%
global H;
global TM;
global D;
global Q;
%
% Compute the function values.
%
fvec=zeros(length(TM),1);
for i=1:length(TM)
  fvec(i)=(Q*exp(-D^2*p(1)/(4*p(2)*TM(i)))/(4*pi*p(2)*TM(i)) - H(i));
end  
