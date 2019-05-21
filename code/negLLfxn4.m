function [ Jtrue ] = negLLfxn4( phat, mudata, vardata, N, t, N0, V0 )
pfxform4 = @(pval)[1 1 1 1].*log(pval); %'forward' parameter transform into Reals
pbxform4 = @(phat)[1 1 1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp

modelfun_mu4 = @(p)mu_fxn(p(1:2), t, N0); % for fitting with sigmas
modelfun_V4= @(p)V_fxn(p(1:2),t,N0, V0);
var_in_mean4 = @(p)(1/N).*(V_fxn(p(1:2), t, N0, V0)) + ((p(3).^2)).*ones(length(t),1);
var_in_var4=@(p)(1/N).*(V4_fxn_ODE(p(1:2),N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p(1:2), t, N0, V0).^2))) + ((p(4).^2).*ones(length(t),1));




J4= @(phat)((0.5*sum((log(2*pi*var_in_mean4(pbxform4(phat)))) + (((yfxform(modelfun_mu4(pbxform4(phat))) - yfxform(mudata))./sqrt(var_in_mean4(pbxform4(phat)))).^2))) +...
 (0.5*sum((log(2*pi*var_in_var4(pbxform4(phat)))) + (((yfxform(modelfun_V4(pbxform4(phat))) - yfxform(vardata))./sqrt(var_in_var4(pbxform4(phat)))).^2))));

objfun_J = @(phat)(J4(phat));
J = objfun_J(phat);
Jtrue = real(J);



end

