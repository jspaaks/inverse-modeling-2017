function [newgen,accept] = metrop_dr(x_new_2,p_new_2,log_p_new_2,x_new,x_old,p_old,alpha12,alpha32,ii,MCMCPar,iR,newgen,N,...
    Sigma,Extra,option);
% Delayed rejection acceptance rule

% Calculate the number of Chains
NrChains = size(p_new_2,1);

% Loop over the chains
for qq = 1:NrChains,
    % Compute l2
    pnew2 = p_new_2(qq,1); pold = p_old(qq,1); logpnew2 = log_p_new_2(qq,1);
    
    % The which option is used
    if option == 1 % Direct probability evaluation
        l2 = (pnew2/pold);
    end;

    if option == 2 | option == 4 | option == 7 | option == 8 % Lnp probability evaluation
        l2 = exp(pnew2 - pold);
    end;

    if option == 3, % SST probability evaluation
        l2 = (pnew2/pold).^(-N.*(1+MCMCPar.Gamma)./2);
    end;
   
    if option == 5, % Special probability evaluation for case study 2
        l2 = exp(-0.5 * (-pnew2 + pold)/Sigma^2); % signs are different because we write -SSR
    end;
    
    % Calculate q1
    q1 = exp(-0.5*(norm((x_new_2(qq,1:MCMCPar.n)-x_new(qq,1:MCMCPar.n))*iR)^2 - norm((x_old(qq,1:MCMCPar.n)-x_new(qq,1:MCMCPar.n))*iR)^2));
    % Calculate alpha            
    alpha13 = l2*q1*(1-alpha32(qq,1))/(1-alpha12(qq,1));

    % Check if prior distribution is being used
    if strcmp(Extra.InitPopulation,'PRIOR');
        % Compute the prior 
        for zz = 1:MCMCPar.n,
            % Define x_prior of proposal 
            x_prior = x_new_2(qq,zz); prior_prop(1,zz) = eval(Extra.prior_dens);
            % Define x_prior of current location chain (old point)
            x_prior = x_old(qq,zz); prior_old(1,zz) = eval(Extra.prior_dens);
        end;

        % Take the product of the various densities
        prior_old = prod(prior_old); prior_prop = max(prod(prior_prop),1e-299); % (to avoid 0/0 --> NaN)
        % Take the ratio
        alpha_pr = prior_prop./prior_old;
        % Now update alpha value with prior
        alpha13 = alpha13.*alpha_pr;  
    end;

    % Generate a random number
    Z = rand;
    
    % Now do MH evaluation    
    if alpha13 > Z
        % accept the new point
        accept(ii(qq),1) = 1;
        % Accept the new point
        newgen(ii(qq),1:MCMCPar.n+2) = [x_new_2(qq,1:MCMCPar.n) pnew2 logpnew2];
    else
        accept(ii(qq),1) = 0;
        % Do nothing
    end;    
end;