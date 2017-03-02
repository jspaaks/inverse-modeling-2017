
% define the number of samples in the population
nPop = 50;

% define the number of dimensions ( =length of the parameter vector)
nDims = 2;

% define which columns in parents and in proposals hold the parameter vector
parCols = 1:nDims;

% preallocate the space needed by the samples and their corresponding
% objective score
proposals = repmat(NaN, [nPop, nDims + 1]);

% iterate over the members of the population
for iPop = 1:nPop

    % draw 4 integer numbers from 1 to nPop
    v = randperm(nPop, 4);
    
    % if iPop happens to be drawn, remove it
    v(v == iPop) = [];
    
    % retrieve selected rows from the parents array; use only the parameter
    % vector columns
    r1 = parents(v(1), parCols);
    r2 = parents(v(2), parCols);
    r3 = parents(v(3), parCols);

    % calculate the two distances dist1 and dist2    
    dist1 = r1 - parents(iPop, parCols);
    dist2 = r3 - r2;

    % calculate the proposal point
    proposals(iPop, parCols) = parents(iPop, parCols) + F * dist1 + K * dist2;
end
