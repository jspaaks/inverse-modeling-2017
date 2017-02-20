function ySim = quadraticmodel(pars,Extra)


a = pars(1);
b = pars(2);
c = pars(3);

xSim = Extra.xObs;

ySim = a*xSim.^2+b*xSim+c;


