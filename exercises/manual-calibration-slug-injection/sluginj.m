function headSim=sluginj(p,TS)
%
% headSim = sluginj(p,TS)
%
% Computes the forward model prediction for the slug test example.
%


% global variables.
global D;
global Q;

% Compute the head:

headSim=repmat(NaN,size(TS));
for i=1:numel(TS)
  headSim(i)=(Q*exp(-D^2*p(1)/(4*p(2)*TS(i)))/(4*pi*p(2)*TS(i)));
end  
