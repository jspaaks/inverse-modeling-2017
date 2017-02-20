function [p,log_p] = CompDensity(x,MCMCPar,Measurement,ModelName,Extra,option)
% This function computes the density of each x value

% Sequential evaluation
for ii=1:size(x,1),

    % Call model to generate simulated data
    evalstr = ['ModPred = ',ModelName,'(x(ii,:),Extra);']; eval(evalstr);

    if option == 1, % Model returns posterior density (mathematical test functions)
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = log(p(ii,1));
    end;

    if option == 2, % Model returns vector of predictions
        % Calculate the error residual
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the log likelihood
        if size(Measurement.Sigma,1) == 1,        
            log_p(ii,1) = N.*log(MCMCPar.Wb./Measurement.Sigma) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma)).^(2/(1+MCMCPar.Gamma))));
        else
            log_p(ii,1) = sum(log(MCMCPar.Wb./Measurement.Sigma)) - MCMCPar.Cb.*(sum((abs(Err./Measurement.Sigma)).^(2/(1+MCMCPar.Gamma))));
        end;
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii];
    end;

    if option == 3, % Model returns vector of predictions
        % Calculate the error residual
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 4, % Model returns log posterior density
        p(ii,1:2) = [ModPred ii]; log_p(ii,1) = p(ii,1);
    end;

    if option == 5, % Similar as option 3, but residuals now weighted with error variance --> homoscedastic measurement error (weighting done in metrop.m)
        % Calculate the error residual
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the sum of squared error
        SSR = sum(abs(Err).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;
    
    if option == 6, % Similar as option 5, but each residual weighted with different error variance --> heteroscedastic measurement error
        % Calculate the error residual
        Err = (Measurement.MeasData(:)-ModPred(:));
        % Derive the weighted sum of squared error
        SSR = sum(abs(Err./Measurement.Sigma).^(2/(1+MCMCPar.Gamma)));
        % And retain in memory
        p(ii,1:2) = [-SSR ii]; log_p(ii,1) = -0.5 * SSR;
    end;

    if option == 7, % Log likelihood with AR-1 model of residuals
        % Calculate the error residual
        Err = (Measurement.MeasData(:)-ModPred(:));
        % First order autoregressive (AR-1) correction of residuals
        rho = x(ii,MCMCPar.n); Err_2 = Err(2:Measurement.N,1) - rho * Err(1:Measurement.N-1,1);
        % Now compute the log-likelihood
        if size(Measurement.Sigma,1) == 1,        
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log(Measurement.Sigma^2 / (1-rho^2)) - (1/2) * (1-rho^2) * (Err(1)/Measurement.Sigma)^2 ...
            - (1/2) * sum((Err_2./Measurement.Sigma).^2);
        else
            log_p(ii,1) = -(Measurement.N/2) * log(2*pi) - (Measurement.N/2) * log(mean(Measurement.Sigma.^2) / (1-rho^2)) - (1/2) * (1-rho^2) * (Err(1)/Measurement.Sigma(1))^2 ...
            - (1/2) * sum((Err_2./Measurement.Sigma(2:Measurement.N)).^2);
        end;        
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii];
    end;

    if option == 8, % Generalized log likelihood (GL)
        % Extract statistical model parameters
        par = Extra.fpar;               % fixed parameters
        par(Extra.idx_vpar) = x(ii,:);  % variable parameters
        par = par';                     % make it a column vector
        statpar = par(end-10:end);
        % Compute the log-likelihood
        log_p(ii,1) = GL('est',statpar,ModPred,Measurement.MeasData);
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii];
    end;
	
	if option == 9, % Limits of acceptability approach (GLUE)
        % Check whether ModPred is in limits of acceptability interval
		idx_1 = ModPred(:) > (Measurement.MeasData-Measurement.Limits);
        idx_2 = ModPred(:) < (Measurement.MeasData+Measurement.Limits);
		% Compute the log_likelihood
		log_p(ii,1) = length(find(idx_1 + idx_2 == 2));
        % And retain in memory
        p(ii,1:2) = [log_p(ii,1) ii];
    end;	
	
end;
