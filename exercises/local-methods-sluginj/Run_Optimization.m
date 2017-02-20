% plot the initial estimate
function [P_GN,P_LM]=Run_Optimization(p_ini)

Exercise1;


% Calculate the objective function as a function of S and T

S=[0:0.002:0.2];
T=[0:0.02:2];

for i=1:101
    for j=1:101
        S1=S(i);
        T1=T(j);
        OF(i,j)=sum(residuals([S1 T1]).^2);
    end
end

% Gauss-Newton method

count=1;
p(count,1:2)=[p_ini];
p(count,3)=sum(residuals(p(count,1:2)).^2);
for count=2:20
J=jac(p(count-1,1:2));
dx=-1*inv(J'*J)*J'*(residuals(p(count-1,:)));
p(count,1:2)=[p(count-1,1:2)+dx'];
p(count,3)=sum(residuals(p(count,1:2)).^2);
end
 
P_GN=p

figure
subplot(221)
pcolor(S,T,OF')
caxis([0 0.5])
colorbar
shading('interp')
xlabel('S')
ylabel('T')
title('Gauss-Newton');
hold on
plot(p(:,1),p(:,2),'m--o');
% axis square

% Uncertainty and correlation in model parameters
DF=6-2;
J=jac(P_GN(end,:));
s2= sum(residuals(P_GN(end,:)).^2)/DF;
unc=s2*inv(J'*J);
S_CL=2.776*sqrt(unc(1,1));
T_CL=2.776*sqrt(unc(2,2));
rho_ST=unc(2,1)./sqrt(unc(1,1)*unc(2,2));

% Levenberg-Marquardt 
clear lambda;
lambda=0.001;
count=1;
clear p;
p(count,1:2)=[p_ini];
p(count,3)=sum(residuals(p(count,:)).^2);
OF_old=sum(residuals(p(count,:)).^2);
for kk=2:50
J=jac(p(count,1:2));
dx=-1*inv(J'*J+lambda*eye(2))*J'*(residuals(p(count,:)));
prop(1:2)=[p(count,1:2)+dx'];
OF_new=sum(residuals(prop).^2);
if (OF_new < OF_old)
    count=count+1;
    p(count,1:3)=[p(count-1,1:2)+dx',OF_new];
    lambda=lambda/2;
    OF_old=OF_new;
   if (lambda <10^(-12))
        lambda=1.0e-12;
   end
else
     count=count+1;
     p(count,1:3)=[p(count-1,:)];
     lambda=lambda*2;
      if (lambda >10^16)
        lambda=10^16;
      end
end
end
P_LM=p

subplot(222)
pcolor(S,T,OF')
caxis([0 0.5])
colorbar
shading('interp')
xlabel('S')
ylabel('T')
title('Levenberg-Marquardt');
hold on
plot(p(:,1),p(:,2),'m--o');
%axis square
