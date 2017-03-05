function objScore = doublegauss(x)

parMu1 = -10;
parSigma1 = 3;
parMu2 = 5;
parSigma2 = 1;


objScore = -1 * (normpdf(x,parMu1,parSigma1) / 2 + ...
                 normpdf(x,parMu2,parSigma2) / 2);
