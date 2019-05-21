function [ Jtrue ] = negLLfxnA( phat, mudata, vardata, N, t, N0, V0 )

pfxform = @(pval)[1 1 1 ].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1 ].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu= @(p)mu_fxnA(p,t, N0); % <n> for params p = b,d, & A
modelfun_V= @(p)V_fxnA(p,t, N0, V0); % variance as a function of parameters
var_in_mean =  @(p)(1/N).*(V_fxnA(p, t, N0, V0)) + (0.01./mudata).^2; % variance in mean
var_in_var = @(p)(1/N).*(V4_fxnA(p,t, N0, V0)-(((N-3)./(N-1)).*(V_fxnA(p, t, N0, V0).^2))) + (0.01./vardata).^2;



% Don't fit on first time point because variance will be 0
J= @(phat)((0.5*sum((log(2*pi*var_in_mean(pbxform(phat)))) + (((yfxform(modelfun_mu(pbxform(phat))) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
 (0.5*sum((log(2*pi*var_in_var(pbxform(phat)))) + (((yfxform(modelfun_V(pbxform(phat))) - yfxform(vardata))./sqrt(var_in_var(pbxform(phat)))).^2))));

objfun_J = @(phat)(J(phat));
J = objfun_J(phat);
Jtrue = real(J);
end