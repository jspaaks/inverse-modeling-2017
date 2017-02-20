% ----------------- DiffeRential Evolution Adaptive Metropolis algorithm ---------------------- %
%                                                                                               %
% DREAM runs multiple different chains simultaneously for global exploration, and automatically %
% tunes the scale and orientation of the proposal distribution using differential evolution.    %
% The algorithm maintains detailed balance and ergodicity and works well and efficient for a    %
% large range of problems, especially in the presence of high-dimensionality and                %
% multimodality. 							                                                    %                     
%                                                                                               %
% DREAM developed by Jasper A. Vrugt and Cajo ter Braak                                         %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%      input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain     %
%      Monte Carlo simulation, Water Resources Research, 44, W00B09, doi:10.1029/2007WR006720,  %
%      2008.                                                                                    %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by differential evolution with         %
%       self-adaptive randomized subspace sampling, International Journal of Nonlinear          %
%       Sciences and Numerical Simulation, 10(3), 271-288, 2009.                                %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson, Equifinality of formal        %
%       (DREAM) and informal (GLUE) Bayesian approaches in hydrologic modeling?, Stochastic     %
%       Environmental Research and Risk Assessment, 1-16, doi:10.1007/s00477-008-0274-y, 2009,  %
%       In Press.                                                                               %   
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution          %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   Ter Braak, C.J.F., and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008.            %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, and J.M. Hyman, Differential evolution adaptive Metropolis   %
%       with snooker update and sampling from past states, SIAM journal on Optimization, 2009.  %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, and J.M. Hyman, Parallel Markov chain Monte Carlo simulation %
%       on distributed computing networks using multi-try Metropolis with sampling from past    %
%       states, SIAM journal on Scientific Computing, 2009.                                     %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
% All rights reserved.                                                                          %
%                                                                                               %
% Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.      %
% Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is     %
% operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.     %
% Government has rights to use, reproduce, and distribute this software.                        %
%                                                                                               %
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR  %
% IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to   %
% produce derivative works, such modified software should be clearly marked, so as not to       %
% confuse it with the version available from LANL.                                              %
%                                                                                               %
% Additionally, redistribution and use in source and binary forms, with or without              %
% modification, are permitted provided that the following conditions are met:                   %
% • Redistributions of source code must retain the above copyright notice, this list of         %
%   conditions and the following disclaimer.                                                    %
% • Redistributions in binary form must reproduce the above copyright notice, this list of      %
%   conditions and the following disclaimer in the documentation and/or other materials         %
%   provided with the distribution.                                                             %
% • Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL %
%   the U.S. Government, nor the names of its contributors may be used to endorse or promote    %
%   products derived from this software without specific prior written permission.              %
%                                                                                               %
% THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND   %
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      %
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS %
% ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF   %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        %
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       %
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                            %
%                                                                                               %
% MATLAB code written by Jasper A. Vrugt, Center for NonLinear Studies (CNLS)                   %
%                                                                                               %
% Written by Jasper A. Vrugt: jasper@uci.edu                                                    %
%                                                                                               %
% Version 0.5: June 2008                                                                        %
% Version 1.0: October 2008       Adaption updated and generalized CR implementation            %
% version 1.1: January 2010       Restart run and new AR-1 likelihood function with model error %
% version 1.2: August 2010        Generalized likelihood function and prior distribution        %
% version 1.3: September 2010     Explicit treatment of prior distribution                      %
% version 1.4: October 2010       Limits of acceptability (GLUE) and few small changes          %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------------- %
%                                                                                               %
%           Note: DE-MC of ter Braak et al. (2006) is a special variant of DREAM                %
%                     It can be executed using the following options                            %
%                                                                                               %
%       MCMCPar.seq = 2 * MCMCPar.n;                                                            %
%       MCMCPar.DEpairs = 1;                                                                    %
%       MCMCPar.nCR = 1;                                                                        %
%       MCMCPar.eps = 0;                                                                        %
%       MCMCPar.outlierTest = 'None';                                                           % 
%       Extra.pCR = 'Fixed'                                                                     %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------------- %
%                                                                                               %
%   After termination of DREAM you can generate a 2D matrix with samples using the command:     %
%                                                                                               %
%       ParSet = GenParSet(Sequences,MCMCPar); if Extra.save_in_memory = 'Yes'                  %
%       otherwise                                                                               %
%       ParSet = GenParSet(Reduced_Seq,MCMCPar); if Extra.reduced_sample_collection = 'Yes'     %
%                                                                                               %
%       output.R_stat   The Gelman Rubin convergence diagnostic                                 %
%       output.AR       The acceptance percentage of candidate points                           %
%       output.CR       Optimized probabilities for individual crossover rates                  %
%       output.outlier  Outlier chain numbers                                                   %
%       hist_logp       Log density of last 50% of samples in each chain                        %
%       X               Final position of chains (population)                                   %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %


% Different test examples
% example 1: n-dimensional banana shaped Gaussian distribution
% example 2: n-dimensional Gaussian distribution
% example 3: n-dimensional multimodal mixture distribution
% example 4: real-world example using hymod model
% example 5: rainfall-runoff model with generalized log-likelihood function
% example 6: HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters
% example 7: HYMOD rainfall - runoff model -- Limits of Acceptability Approach (GLUE)

% --------------------------------------------------------------------------------------------- %
%                                                                                               %
%           Note: PRIOR INFORMATION IN DREAM CAN BE SPECIFIED AS FOLLOWS                        %
%                                                                                               %
%       Set Extra.InitPopulation = 'PRIOR'. Then for each parameter specify the prior           %
%       distribution using MATLAB language.                                                     %
%                                                                                               %
%       Example: Four parameters: (1) weibull, (2) lognormal, (3) normal distribution,          %
%       and (4) uniform                                                                         % 
%                                                                                               %
%       Extra.prior = {'wblrnd(9,3)','lognrnd(0,2)','normrnd(-2,3)','unifrnd(0,10)'};           %
%                                                                                               %
%       Weibull:   scale = 9; shape = 3;                                                        %
%       Lognormal: mu = 0; sigma = 2;                                                           %
%       Normal:    mu = -2; sigma = 3;                                                          %
%       Uniform:   lower bound = 0; upper bound = 10;                                           %
%                                                                                               %
%       To calculate the Metropolis ratio, DREAM will assume that the respective densities      %
%       of the various prior distributions end with "pdf" rather than "rnd"                     %                                %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% clear
% close all
% clc

addpath('.\..\..\dream')

% Which example to run?
example = 3

if example == 1, % n-dimensional banana shaped Gaussian distribution

    MCMCPar.n = 10;                         % Dimension of the problem (Nr. parameters to be optimized in the model)
    MCMCPar.seq = 10;                       % Number of Markov Chains / sequences
    MCMCPar.DEpairs = 3;                    % Number of chain pairs to generate candidate points
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.ndraw = 25000;                 % Maximum number of function evaluations
    MCMCPar.steps = 10;                     % Number of steps
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 2;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % Define the specific properties of the banana function
    Extra.mu   = [zeros(1,MCMCPar.n)];                      % Center of the banana function
    Extra.cmat = eye(MCMCPar.n); Extra.cmat(1,1) = 100;     % Target covariance
    Extra.imat = inv(Extra.cmat);                           % Inverse of target covariance
    Extra.bpar = [0.1];                                     % "bananity" of the target, see bananafun.m

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = Extra.mu;                                   % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n) * 5;                        % Initial covariance
    % Save all information in memory or not?
    Extra.save_in_memory = 'Yes'; % was 'No'

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];

    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'Banshp';
    % Define likelihood function
    option = 4;

    % specify the visualization routine
    Extra.visScriptName = 'visDream';

end;

if example == 2,    % n-dimensional Gaussian distribution

    MCMCPar.n = 10; %was 100                        % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 1000000;                % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 3;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 10;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-5 * ones(1,MCMCPar.n)]; ParRange.maxn = [15 * ones(1,MCMCPar.n)];
    %ParRange.minn = [9.9 * ones(1,MCMCPar.n)]; ParRange.maxn = [10 * ones(1,MCMCPar.n)];

    % Define the boundary handling
    Extra.BoundHandling = 'None';

    % ---------------------- Define covariance matrix ---------------------
    % Construct the d x d covariance matrix
    A = 0.5*eye(MCMCPar.n) + 0.5*ones(MCMCPar.n);
    % Rescale to variance-covariance matrix of interest
    for i=1:MCMCPar.n
        for j=1:MCMCPar.n
            C(i,j) = A(i,j)*sqrt(i*j);
        end
    end
    % Set to Extra
    Extra.qcov = C; Extra.muX = zeros(1,MCMCPar.n); Extra.invC = inv(C);
    % ---------------------------------------------------------------------

    Extra.save_in_memory = 'Yes';  % was 'No'
    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define center of target
    Extra.mu = zeros(1,MCMCPar.n);          
    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Define modelName
    ModelName = 'normalfunc';
    % Define likelihood function
    option = 4;

end;

if example == 3,    % n-dimensional multimodal mixture distribution

    MCMCPar.n = 10;                         % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 1000000;                % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 3;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 5e-2;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'Yes';% Thinned sample collection?
    Extra.T = 10;                           % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'COV_BASED';
    % Provide information to do alternative sampling
    Extra.muX = zeros(1,MCMCPar.n);         % Provide mean of initial sample
    Extra.qcov = eye(MCMCPar.n);            % Initial covariance

    Extra.Lam = eye(MCMCPar.n);             % covariance
    Extra.mu1 = -5 * ones(1,MCMCPar.n);     % center point of first density
    Extra.mu2 =  5 * ones(1,MCMCPar.n);     % center point of second density
    Extra.sigma = eye(MCMCPar.n);           % define the std of the target

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [-Inf * ones(1,MCMCPar.n)]; ParRange.maxn = [Inf * ones(1,MCMCPar.n)];
    % Define the boundary handling
    Extra.BoundHandling = 'None';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';  % was 'No'

    % Define data structures for use in computation of posterior density
    Measurement.MeasData = []; Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'mixturemodel';
    % Define likelihood function
    option = 1;

end;

if example == 4,    % HYMOD rainfall - runoff model

    MCMCPar.n = 5;                          % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = MCMCPar.n*1000;         % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 1;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 2e-1;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10];
    ParRange.maxn = [500 2.00 0.99 0.10 0.99];
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the Leaf River data
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795; 

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.PET = bound(1:Extra.MaxT,5); Extra.Precip = sum(bound(1:Extra.MaxT,6:9),2);

    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'hymod';
    % Define likelihood function -- Sum of Squared Error
    option = 2; Measurement.Sigma = ones(Measurement.N,1);

    % specify the visualization routine
    Extra.visScriptName = 'visDream';
    % Specify labels for the parameter names
    Extra.ParNames.labels = {'c_{max}','b_{exp}','\alpha','R_s','R_q'};
    Extra.ParNames.interpreter = 'tex';


end;

if example == 5,    % Rainfall-runoff model with generalized log-likelihood

    % ---------------------------- Check the following 2 papers ------------------------------

    % G. Schoups, J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption of 
    %     accuracy and efficiency of Markov Chain Monte Carlo simulation by inaccurate 
    %     numerical implementation of conceptual hydrologic models, Water Resources
    %     Research, 46, W10530, doi:10.1029/2009WR008648, In Press. 

    % G. Schoups, and J.A. Vrugt (2010), A formal likelihood function for parameter and 
    %     predictive inference of hydrologic models with correlated, heteroscedastic and 
    %     non-Gaussian errors, Water Resources Research, 46, W10531, doi:10.1029/2009WR008933. 

    % ----------------------------------------------------------------------------------------

    MCMCPar.n = 11;                         % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 30000;                  % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 1;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 2e-1;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';
    % Give the parameter ranges (minimum and maximum values)
    %parno:       1     2     3     4     5     6     7      8    9     10    11   12   13   14   15   16   17   18   19   20  21
    %parname:     fA    Imax  Smax  Qsmax alE   alS   alF    Pf   Kfast Kslow std0 std1 beta xi   mu1  phi1 phi2 phi3 phi4 K   lambda
    Extra.fpar = [1     0     100   10    100   1e-6  1e-6   0    2     70    0.1  0    0    1    0    0    0    0    0    0   1];
    parmin =     [1     0     10    0     1e-6  1e-6 -10     0    0     0     0    0   -1    0.1  0    0    0    0    0    0   0.1 ];
    parmax =     [1     10    1000  100   100   1e-6  10     0    10    150   1    1    1    10   100  1    1    1    1    1   1];
    Extra.idx_vpar = [2 3 4 5 7 9 10 11 12 13 16];    
    ParRange.minn = parmin(Extra.idx_vpar); ParRange.maxn = parmax(Extra.idx_vpar);
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the French Broad data
    daily_data = load('03451500.dly');

    % First two years are warm-up
    Extra.idx = [731:size(daily_data,1)]';

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.Precip = daily_data(:,4);
    Extra.Ep     = daily_data(:,5);

    % Define the measured streamflow data
    Measurement.MeasData = daily_data(Extra.idx,6); 
    Measurement.Sigma = []; 
    Measurement.N = size(Measurement.MeasData,1);
    % Define modelName
    ModelName = 'hmodel';
    % Use generalized likelihood function
    option = 8; 
end;

if example == 6,	% HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters

    % ---------------------------- Check the following paper ------------------------------
    %   B. Scharnagl, J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse
    %	modeling of soil water dynamics at the field scale: using prior information
    %	on soil hydraulic properties, Hydrology and Earth System Sciences.
    % -------------------------------------------------------------------------------------

    MCMCPar.n = 7;                          % Dimension of the problem
    MCMCPar.seq = 10;                      % Number of Markov Chains / sequences
    MCMCPar.ndraw = 5000;                   % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 1;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 2e-1;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';

    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Define feasible parameter space (minimum and maximum values)
    %				1		2		3				4			5			6		7     
    %				[thetar	thetas	log10(alpha)	log10(n)	log10(Ks)	L		hLB
    ParRange.minn =	[0.0430	0.4090	-2.5528			0.1790		-2.2366		-5.4900	-250];
    ParRange.maxn =	[0.0910 0.4810	-2.0706			0.2670		-0.0800		6.2700	-50];

    % Store working directory and subdirectory containing the files needed to run this example
    Extra.workdir = pwd;
    Extra.subdir = [pwd '\example_' num2str(example)];

    % Add subdirectory to search path
    addpath(Extra.subdir)

    % Provide observational data and data needed to modify the initial and boundary conditions
    [Measurement,Extra] = ProvideData(Extra);

    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';

    % Define model name
    ModelName = 'HYDRUS';

    % Define option (model computes log-likelihood)
    option = 4;

    % Indicate the use prior information
    Extra.InitPopulation = 'PRIOR';

    % Specify the prior distributions for the various parameters
    Extra.prior = {'normrnd(0.0670,0.0060)',...
                   'normrnd(0.4450,0.0090)',...
                   'normrnd(-2.310,0.0600)',...
                   'normrnd(0.2230,0.0110)',...
                   'normrnd(-1.160,0.2700)',...
                   'normrnd(0.3900,1.4700)',...
                   'unifrnd(-250,-50)'};

    Extra.ParNames.labels = {'\theta_{r}','\theta_{s}','log10(\alpha)','log10(n)','log10(K_{s})','L','hLB'};
    Extra.ParNames.interpreter = 'tex';

end;


if example == 7,    % HYMOD rainfall - runoff model -- Limits of Acceptability Approach (GLUE)

    MCMCPar.n = 5;                          % Dimension of the problem
    MCMCPar.seq = MCMCPar.n;                % Number of Markov Chains / sequences
    MCMCPar.ndraw = 15000;                   % Maximum number of function evaluations
    MCMCPar.nCR = 3;                        % Crossover values used to generate proposals (geometric series)
    MCMCPar.Gamma = 0;                      % Kurtosis parameter Bayesian Inference Scheme
    MCMCPar.DEpairs = 1;                    % Number of DEpairs
    MCMCPar.steps = 10;                     % Number of steps in sem
    MCMCPar.eps = 2e-1;                     % Random error for ergodicity
    MCMCPar.outlierTest = 'IQR_test';       % What kind of test to detect outlier chains?
    MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

    % -----------------------------------------------------------------------------------------------------------------------
    Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
    % -----------------------------------------------------------------------------------------------------------------------

    % --------------------------------------- Added for reduced sample storage ----------------------------------------------
    Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
    Extra.T = 1000;                         % Every Tth sample is collected
    % -----------------------------------------------------------------------------------------------------------------------

    % What type of initial sampling
    Extra.InitPopulation = 'LHS_BASED';

    % Give the parameter ranges (minimum and maximum values)
    ParRange.minn = [1.0 0.10 0.10 0.00 0.10]; 
    ParRange.maxn = [500 2.00 0.99 0.10 0.99];
    % Define the boundary handling
    Extra.BoundHandling = 'Reflect';
    % Save in memory or not
    Extra.save_in_memory = 'Yes';

    % Load the Leaf River data
    load bound.txt;

    % Then read the boundary conditions -- only do two years
    Extra.MaxT = 795; 

    % Define the PET, Measured Streamflow and Precipitation.
    Extra.PET = bound(1:Extra.MaxT,5); Extra.Precip = sum(bound(1:Extra.MaxT,6:9),2);

    % Define the measured streamflow data
    Measurement.MeasData = bound(65:Extra.MaxT,4); Measurement.Sigma = []; Measurement.N = size(Measurement.MeasData,1);

    % Define modelName
    ModelName = 'hymod';
    % Define likelihood function -- Limits of acceptability approach (GLUE)
    option = 9; 
    % Define the limits
    Measurement.Limits = 0.1 * Measurement.MeasData;

    % specify the visualization routine (empty string for no visualization)
    Extra.visScriptName = 'visDream';
    % Specify labels for the parameter names
    Extra.ParNames.labels = {'c_{max}','b_{exp}','\alpha','R_s','R_q'};
    Extra.ParNames.interpreter = 'tex';

end;

% Scale of Delayed Rejection if used in algorithm
Extra.DR = 'No'; Extra.DRscale = 10; Restart = 'No';

% Run the MCMC algorithm
[Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);


figure(1)
print(['example-',num2str(example),'-objective-scores.eps'],'-r300','-loose','-depsc2')
print(['example-',num2str(example),'-objective-scores.png'],'-r300','-dpng')
figure(2)
print(['example-',num2str(example),'-sequences.eps'],'-r300','-loose','-depsc2')
print(['example-',num2str(example),'-sequences.png'],'-r300','-dpng')

