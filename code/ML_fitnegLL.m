function [negLLfit,pbest] = ML_fitnegLL(theta, tsamp, Ninit, N,mudatavec, vardatavec, modelcode, profindex, currparams)

switch modelcode
    case 1
        
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
    
            if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_donly = @(d) objfun_J(pfxform([currb,d]));
            lb = [0 ];
            ub = [ Inf ];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_donly, theta(2), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bonly = @(b) objfun_J(pfxform([b,currd]));
            lb = [0 ];
            ub = [ Inf ];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bonly, theta(1), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
    case 2
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));


        % PERFORM FITTING OF ALL DATA TO SINGLE EXPONENTIAL MODEL
        pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
        pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
        yfxform = @(y)log(y); % 'forward' transform for data and model output
        ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
        
        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

            if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 0];
            ub = [ Inf  Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAonly, theta(2:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAonly = @(p) objfun_J(pfxform([p(1),currd, p(2)]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAonly, [theta(1), theta(3)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdonly = @(p) objfun_J(pfxform([p,currA]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdonly, theta(1:2), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end


case 3  % strong Allee model Allee death
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));
        
       
        pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
        pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
        yfxform = @(y)log(y); % 'forward' transform for data and model output
        ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
            sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));

            if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 0];
            ub = [ Inf  Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAonly, theta(2:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAonly = @(p) objfun_J(pfxform([p(1),currd, p(2)]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAonly, [theta(1), theta(3)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdonly = @(p) objfun_J(pfxform([p,currA]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdonly, theta(1:2), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end

        

case 4  % strong Allee model Allee on birth & death
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));
        
       
        pfxform = @(pval)[1 1 1].*log(pval); %'forward' parameter transform into Reals
        pbxform = @(phat)[1 1 1].*exp(phat); %'backward' parameter transform into model space
        yfxform = @(y)log(y); % 'forward' transform for data and model output
        ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output


        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
            sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));


            if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 0];
            ub = [ Inf  Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAonly, theta(2:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAonly = @(p) objfun_J(pfxform([p(1),currd, p(2)]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAonly, [theta(1), theta(3)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdonly = @(p) objfun_J(pfxform([p,currA]));
            lb = [0 0];
            ub = [ Inf Inf];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdonly, theta(1:2), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end



case 5 % weak Allee model Allee on birth
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));
        
        lb = [0 0 -10 -10];
        ub = [ 2 2 10 10];
       
        pfxform = @(pval)[1 1 1 1].*(pval-lb)./(ub-lb); %'forward' parameter transform (normalized)
        pbxform = @(phat) [1 1 1 1].*((ub-lb).*(phat))+lb;
        %yfxform = @(y)log(y); % 'forward' transform for data and model output
        %ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
        


        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
            sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
            if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAtonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 -10 0 ];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAtonly, theta(2:4), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAtonly = @(p) objfun_J(pfxform([p(1),currd, p(2:3)]));
            lb = [0 -10 0];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAtonly, [theta(1), theta(3:4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdtonly = @(p) objfun_J(pfxform([p(1:2),currA, p(3)]));
            lb = [0 0 0 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdtonly, [theta(1:2), theta(4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==4
            currtau = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdAonly = @(p) objfun_J(pfxform([p(1:3),currtau]));
            lb = [0 0 -10 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdAonly, theta(1:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end

        


case 6 % weak Allee model Allee on death
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));
        
        lb = [0 0 -10 0];
        ub = [ 2 2 10 10];
       
        pfxform = @(pval)[1 1 1 1].*(pval-lb)./(ub-lb); %'forward' parameter transform (normalized)
        pbxform = @(phat) [1 1 1 1].*((ub-lb).*(phat))+lb;
        %yfxform = @(y)log(y); % 'forward' transform for data and model output
        %ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
        


        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
            sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
        if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAtonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 -10 0 ];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAtonly, theta(2:4), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAtonly = @(p) objfun_J(pfxform([p(1),currd, p(2:3)]));
            lb = [0 -10 0];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAtonly, [theta(1), theta(3:4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdtonly = @(p) objfun_J(pfxform([p(1:2),currA, p(3)]));
            lb = [0 0 0 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdtonly, [theta(1:2), theta(4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==4
            currtau = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdAonly = @(p) objfun_J(pfxform([p(1:3),currtau]));
            lb = [0 0 -10 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdAonly, theta(1:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end

    


case 7 % weak Allee model Allee on birth & death
        modelfun_mu =@(p)gen_model_mu(p,tsamp, Ninit, modelcode);
        modelfun_V = @(p)gen_model_var(p, tsamp, Ninit, modelcode);
        modelfun_V4= @(p)gen_model_v4(p, tsamp, Ninit, modelcode);
        var_in_mean =  @(p)(1/N).*(modelfun_V(p)); % vertical
        var_in_var = @(p)(1/N).*(modelfun_V4(p)-(((N-3)./(N-1)).*(modelfun_V(p).^2)));
        
        lb = [0 0 -10 0];
        ub = [ 2 2 10 10];
       
        pfxform = @(pval)[1 1 1 1].*(pval-lb)./(ub-lb); %'forward' parameter transform (normalized)
        pbxform = @(phat) [1 1 1 1].*((ub-lb).*(phat))+lb;
        %yfxform = @(y)log(y); % 'forward' transform for data and model output
        %ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output
        


        J = @(phat) (sum(((mudatavec-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
            sum(((vardatavec-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
        
        if profindex==1
            currb = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_dAtonly = @(p) objfun_J(pfxform([currb,p]));
            lb = [0 -10 0 ];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_dAtonly, theta(2:4), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==2
            currd = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bAtonly = @(p) objfun_J(pfxform([p(1),currd, p(2:3)]));
            lb = [0 -10 0];
            ub = [ 2 10 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bAtonly, [theta(1), theta(3:4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==3
            currA = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdtonly = @(p) objfun_J(pfxform([p(1:2),currA, p(3)]));
            lb = [0 0 0 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdtonly, [theta(1:2), theta(4)], options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end
            
            if profindex==4
            currtau = currparams; 
            objfun_J = @(phat)(J(phat));
            objfun_bdAonly = @(p) objfun_J(pfxform([p(1:3),currtau]));
            lb = [0 0 -10 ];
            ub = [ 2 2 10];
            options = optimset('TolX', 1e-8, 'MaxFunEvals', 1000);
            [phatbest_J,fval,exitflag] = fminsearch(objfun_bdAonly, theta(1:3), options);% need a better way to search parameter space,i.e. montecarlo
            negLLfit = fval;
            pbest = exp(phatbest_J);
            end

    
end
