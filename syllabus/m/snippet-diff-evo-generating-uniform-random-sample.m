% define the number of samples in the population
nPop = 50;

% define the number of dimensions ( =length of the parameter vector)
nDims = 2;

% set the lower limits for every dimension of the parameter vector, repeat
% for the number of members in the population
parSpaceLowerBounds = repmat([-100, -100], [nPop,1]);

% set the upper limits for every dimension of the parameter vector, repeat
% for the number of members in the population
parSpaceUpperBounds = repmat([100, 100], [nPop, 1]);

% calculate the range between the lower and upper boundaries on the
% parameter space for each dimension 
parSpaceRange = parSpaceUpperBounds - parSpaceLowerBounds;

% perform nPop * nDims random draws from a uniform distribution
uniformRandomDraw = rand(nPop,nDims);

% fill the parents array with uniform draws within the parameter space
parents(1:nPop,1:nDims) = parSpaceLowerBounds + uniformRandomDraw .* parSpaceRange;
