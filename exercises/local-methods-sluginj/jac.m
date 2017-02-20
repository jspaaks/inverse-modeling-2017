%
%
%
function J=jac(p)
%
% global variables.
%
global H;
global TM;
global D;
global Q;
%
%
%
n=length(TM);
J=zeros(n,2);
for i=1:n
  J(i,1)=(-Q*D^2*exp(-D^2*p(1)/(4*p(2)*TM(i)))/...
         (16*pi*p(2)^2*TM(i)^2));
  J(i,2)=(Q/(4*pi*p(2)^2*TM(i)))*...
          ((D^2*p(1))/(4*p(2)*TM(i))-1)*...
          exp(-D^2*p(1)/(4*p(2)*TM(i)));
end
