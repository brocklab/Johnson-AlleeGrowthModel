function [ Jtrue ] = negLLfxn5( phat, mudata, vardata, N, t, N0, V0 )

pfxform5 = @(pval)[1 1 1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform5 = @(phat)[1 1 1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu5= @(p)mu_fxnA(p(1:3),t, N0); % <n> for params p = b,d, & A
modelfun_V5= @(p)V_fxnA(p(1:3),t, N0, V0); % variance as a function of parameters
var_in_mean5 =  @(p)(1/N).*(V_fxnA(p(1:3), t, N0, V0)) + (p(4).^2).*ones(length(t),1); % variance in mean
var_in_var5 = @(p)(1/N).*(V4_fxnA(p(1:3),t, N0, V0)-(((N-3)./(N-1)).*(V_fxnA(p(1:3), t, N0, V0).^2))) + ((p(5).^2).*ones(length(t),1));



% Don't fit on first time point because variance will be 0
J5= @(phat)((0.5*sum((log(2*pi*var_in_mean5(pbxform5(phat)))) + (((yfxform(modelfun_mu5(pbxform5(phat))) - yfxform(mudata))./sqrt(var_in_mean5(pbxform5(phat)))).^2))) +...
 (0.5*sum((log(2*pi*var_in_var5(pbxform5(phat)))) + (((yfxform(modelfun_V5(pbxform5(phat))) - yfxform(vardata))./sqrt(var_in_var5(pbxform5(phat)))).^2))));

objfun_J5 = @(phat)(J5(phat));
J = objfun_J5(phat);
Jtrue = real(J);
end