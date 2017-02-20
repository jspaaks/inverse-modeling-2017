function ySim = linearmodel(pars,Extra)


a = pars(1);
b = pars(2);

xSim = Extra.xObs;

ySim = a*xSim+b;


