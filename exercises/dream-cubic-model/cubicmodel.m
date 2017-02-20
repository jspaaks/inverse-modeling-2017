function ySim = cubicmodel(pars,Extra)


a = pars(1);
b = pars(2);
c = pars(3);
d = pars(4);

xSim = Extra.xObs;

ySim = a*xSim.^3+b*xSim.^2+c*xSim+d;


