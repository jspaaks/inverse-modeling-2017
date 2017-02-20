%
% fvec=model(p,TM)
%
% Computes the forward model prediction for the slug test example.
%
function fvec=model(p,TM)
%
% global variables.
%
global H;
global D;
global Q;
%
% Compute the function values.
%
fvec=zeros(length(TM),1);
for i=1:length(TM)
  fvec(i)=(Q*exp(-D^2*p(1)/(4*p(2)*TM(i)))/(4*pi*p(2)*TM(i)));
end  
