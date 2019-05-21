modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));


        % PERFORM FITTING OF ALL DATA TO SINGLE EXPONENTIAL MODEL
        pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
        pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
        yfxform = @(y)log(y); % 'forward' transform for data and model output
        ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
        
        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

    objfun_J = @(phat)(J(phat));
    lb = [0 0 ];
    ub = [ Inf Inf ];
    options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
    [phatbest_J,fval,exitflag] = fminsearch(objfun_J, pfxform(theta), options);% need a better way to search parameter space,i.e. montecarlo
    negLLguess= objfun_J(pfxform(theta));
    negLLfit= objfun_J(phatbest_J);
    params_best_J= pbxform(phatbest_J);
    pbest = params_best_J;

    k = 2;
    n = length(mudatavec)*length(Ninit);