function [Sequences,Reduced_Seq,X,output,hist_logp] = dream(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart)
% ----------------- DiffeRential Evolution Adaptive Metropolis algorithm -----------------------%
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
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Check whether to do a restart or not
if strcmp(Restart,'No');

    % Set random number generator
    opts.Seed = 'sum(100*clock)  % evaluated if it is a string';
    % Then generate new seed
    if ischar(opts.Seed)
        randn('state', eval(opts.Seed));     % random number generator state
    else
        randn('state', opts.Seed);
    end

    % Initialize variables used by DREAM
    [MCMCPar,pCR,CR,lCR,Iter,counter,teller,new_teller,hist_logp,Sequences,Table_JumpRate,Reduced_Seq,iloc,output] = InitVariables(MCMCPar,...
        Extra);

    % Step 1: Sample MCMCPar.seq (N) points in the parameter space
    if strcmp(Extra.InitPopulation,'LHS_BASED'),
        % Latin hypercube sampling when indicated
        [x] = LHSU(ParRange.minn,ParRange.maxn,MCMCPar.seq);
    elseif strcmp(Extra.InitPopulation,'COV_BASED');
        % Do alternative sampling method (for some examples)
        [x] = repmat(Extra.muX,MCMCPar.seq,1) + randn(MCMCPar.seq,MCMCPar.n) * chol(Extra.qcov);
    elseif strcmp(Extra.InitPopulation,'PRIOR');
        % Create the initial position of each chain by drawing each parameter individually from the prior
        for qq = 1:MCMCPar.n,
            for zz = 1:MCMCPar.seq,
                x(zz,qq) = eval(char(Extra.prior(qq)));
            end;
        end;
    end;

    % Do boundary handling -- what to do when points fall outside bound
    if strcmp(Extra.BoundHandling,'Reflect');
        [x] = ReflectBounds(x,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'Bound');
        [x] = SetToBounds(x,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'Fold');
        [x] = FoldBounds(x,ParRange);
    end;
    if strcmp(Extra.BoundHandling,'None');
        % Do nothing
    end;

    % Step 2: Calculate posterior density associated with each value in x
    [p,log_p] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra,option);

    % Save the initial population, density and log density in one matrix X
    X = [x p(:,1) log_p];

    % Then initialize the sequences
    if strcmp(Extra.save_in_memory,'Yes');
        [Sequences] = InitSequences(X,Sequences,MCMCPar);
    else
        [Sequences] = InitSequences(X,[],MCMCPar);
    end;

    if strcmp(Extra.reduced_sample_collection,'Yes');
        % Reduced sample collection
        iloc_2 = 0;
    else
        % Do nothing
    end;

    % Save N_CR in memory and initialize delta_tot
    output.CR(1,1:size(pCR,2)+1) = [Iter pCR]; delta_tot = zeros(1,MCMCPar.nCR);

    % Save history log density of individual chains
    hist_logp(1,1:MCMCPar.seq+1) = [Iter reshape(X(:,MCMCPar.n+2),MCMCPar.seq,1)'];

    % Compute the R-statistic
    [output.R_stat(1,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(1:iloc,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];

else % If a restart run is being done: just load the output from the previous ongoing trial

    % Load all information from previously terminated run -- double the number of draws
    load everything.mat; MCMCPar.ndraw = 2 * MCMCPar.ndraw;

end;

% Now start iteration ...
while (Iter < MCMCPar.ndraw),

    % Loop a number of times
    for gen_number = 1:MCMCPar.steps,

        % Initialize DR properties
        accept2 = 0; ItExtra = 0; new_teller = new_teller + 1;

        % Define the current locations and associated posterior densities
        [xold,p_xold,log_p_xold] = GetLocation(X,MCMCPar); ItExtra = 0;

        % Now generate candidate in each sequence using current point and members of X
        [xnew,CR(:,gen_number)] = offde(xold,X,CR(:,gen_number),MCMCPar,Table_JumpRate,ParRange,Extra.BoundHandling,[],'No');

        % Now compute the likelihood of the new points
        [p_xnew,log_p_xnew] = CompDensity(xnew,MCMCPar,Measurement,ModelName,Extra,option); p_xnew = p_xnew(:,1);

        % Now apply the acceptance/rejectance rule for the chain itself
        [newgen,alpha12,accept] = metrop(xnew,p_xnew,log_p_xnew,xold,p_xold,log_p_xold,Measurement.N,Measurement.Sigma,MCMCPar,...
            Extra,option);

        % Check whether we do delayed rejection or not
        if strcmp(Extra.DR,'Yes'),

            % Compute the Cholesky Decomposition of X
            R = (2.38/sqrt(MCMCPar.n)) * chol(cov(xold(1:end,1:MCMCPar.n)) + 1e-5*eye(MCMCPar.n));

            % Compute inverse of R
            iR = inv(R);

            % Scale the adjusted proposal distribution for DR
            R2 = R./Extra.DRscale;

            % Find which chains to actually apply delayed rejection to
            [ii] = find(accept==0); Nii = size(ii,1);

            % Only do DR to selected chains
            if Nii > 0,

                % Now do a delayed rejection step for each chain (but only ii is used later on)
                [xnew2,dummy_2] = offde(xold,[],[],MCMCPar,[],ParRange,Extra.BoundHandling,R2,'Yes');

                % Now compute the likelihood of the new points
                [p_xnew2(ii,1:2),log_p_xnew2(ii,1)] = CompDensity(xnew2(ii,1:MCMCPar.n),MCMCPar,Measurement,ModelName,Extra,option); p_xnew2 = p_xnew2(:,1);

                % Now compute alpha32 (note we turned p_new_2 and p_new around!)
                [dummy_1,alpha32(ii,1),dummy_2] = metrop(xnew(ii,:),p_xnew(ii,1),log_p_xnew(ii,1),xnew2(ii,:),p_xnew2(ii,1),...
                    log_p_xnew2(ii,1),Measurement.N,Measurement.Sigma,MCMCPar,Extra,option);

                % And do evaluation rule
                [newgen,accept2] = metrop_dr(xnew2(ii,1:MCMCPar.n),p_xnew2(ii,1),log_p_xnew2(ii,1),xnew(ii,1:MCMCPar.n),xold(ii,1:MCMCPar.n),...
                    p_xold(ii,1),alpha12(ii,1),alpha32(ii,1),ii,MCMCPar,iR,newgen,Measurement.N,Measurement.Sigma,Extra,option);

                % Update ItExtra
                ItExtra = ItExtra + Nii;

            end;

        end;

        % Define the location in the sequence
        if strcmp(Extra.save_in_memory,'Yes');
            % Define idx based on iloc
            iloc = iloc + 1;
        else
            % Do nothing -- remain at first location
        end

        % Now update the locations of the Sequences with the current locations
        Sequences(iloc,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(newgen',1,MCMCPar.n+2,MCMCPar.seq);

        % Check reduced sample collection
        if strcmp(Extra.reduced_sample_collection,'Yes');
            if (new_teller == Extra.T),
                % Update iloc_2 and new_teller
                iloc_2 = iloc_2 + 1; new_teller = 0;
                % Reduced sample collection
                Reduced_Seq(iloc_2,1:MCMCPar.n+2,1:MCMCPar.seq) = reshape(newgen',1,MCMCPar.n+2,MCMCPar.seq);
            else
                % Do nothing -- not yet at Extra.T
            end;
        end;

        % And update X using current members of Sequences
        X = newgen; clear newgen;

        if strcmp(Extra.pCR,'Update');
            % Calculate the standard deviation of each dimension of X
            r = repmat(std(X(1:MCMCPar.seq,1:MCMCPar.n)),MCMCPar.seq,1);
            % Compute the Euclidean distance between new X and old X
            delta_normX = sum(((xold(1:end,1:MCMCPar.n) - X(1:end,1:MCMCPar.n))./r).^2,2);
            % Use this information to update sum_p2 to update N_CR
            [delta_tot] = CalcDelta(MCMCPar,delta_tot,delta_normX,CR(1:MCMCPar.seq,gen_number));
        end;

        % Update hist_logp
        hist_logp(counter,1:MCMCPar.seq+1) = [Iter X(1:MCMCPar.seq,MCMCPar.n+2)'];

        % Save some important output -- Acceptance Rate
        output.AR(counter,1:2) = [Iter 100 * (sum(accept) + sum(accept2))/(MCMCPar.seq + ItExtra)];
        
        
        if exist('Extra','var')
            if ~isfield(Extra,'visScriptName')
                Extra.visScriptName = 'visDream';
            end
            eval(Extra.visScriptName)
        end
        % Update Iteration and counter
        Iter = Iter + MCMCPar.seq + ItExtra; counter = counter + 1;

    end;

    % ---------------------------------------------------------------------

    % Store Important Diagnostic information -- Probability of individual crossover values
    output.CR(teller,1:size(pCR,2)+1) = [Iter pCR];

    % Do this to get rounded iteration numbers
    if (teller == 2), MCMCPar.steps = MCMCPar.steps + 1; end;

    % Check whether to update individual pCR values
    if (Iter <= 0.1 * MCMCPar.ndraw);
        if strcmp(Extra.pCR,'Update');
            % Update pCR values
            [pCR,lCR] = AdaptpCR(MCMCPar,CR,delta_tot,lCR);
        end;
    else
        % See whether there are any outlier chains, and remove them to current best value of X
        [X,Sequences(iloc,:,:),hist_logp(1:counter-1,2:MCMCPar.seq+1),output.outlier] = RemOutlierChains(X,Sequences(iloc,:,:),...
            hist_logp(1:counter-1,2:MCMCPar.seq+1),Iter,output.outlier,MCMCPar);
    end;

    if strcmp(Extra.pCR,'Update');
        % Generate CR values based on current pCR values
        [CR] = GenCR(MCMCPar,pCR); pCR
    else
        CR = pCR * ones(MCMCPar.seq,MCMCPar.steps);
    end;

    % Calculate Gelman and Rubin convergence diagnostic
    if strcmp(Extra.save_in_memory,'Yes');
        start_loc = max(1,floor(0.5*iloc)); end_loc = iloc;
        % Compute the R-statistic using 50% burn-in from Sequences
        [output.R_stat(teller,1:MCMCPar.n+1)] = [Iter Gelman(Sequences(start_loc:end_loc,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];
    else
        start_loc = max(1,floor(0.5*iloc_2)); end_loc = iloc_2;
        % Compute the R-statistic using 50% burn-in from Reduced_Seq
        [output.R_stat(teller,1:MCMCPar.n+1)] = [Iter Gelman(Reduced_Seq(start_loc:end_loc,1:MCMCPar.n,1:MCMCPar.seq),MCMCPar)];
    end;

    % Update the teller
    teller = teller + 1;

end;

% Postprocess output from DREAM before returning arguments
[i] = find(sum(Sequences(:,1:end,1),2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; Sequences = Sequences(1:i,:,:);
end;

[i] = find(sum(Reduced_Seq(:,1:end),2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; Reduced_Seq = Reduced_Seq(1:i,:,:);
end;

% Postprocess output from DREAM before returning arguments
[i] = find(sum(output.R_stat(:,1:end,1),2)==0);
% Remove extra rows from Sequences
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; output.R_stat = output.R_stat(1:i,:,:);
end;

[i] = find(sum(output.AR(:,1:end),2)==0);
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; output.AR = output.AR(1:i,:);
end;

[i] = find(sum(output.CR(:,1:end),2)==0);
if isempty(i) == 0,
    % Derive i
    i = i(1) - 1; output.CR = output.CR(1:i,:);
end;



