function [ Jtrue ] = garchlike( phat, mudata, vardata, N, t, N0, V0 )
pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

% Two observations to fit on, mu_data and var_data
% both use tsamp
modelfun_mu = @(p)mu_fxn(p, t, N0); % <n> for params (p) = b & d, vertical
modelfun_n2= @(p)second_mom_fxn(p, t, N0); % Var for params (p) = b & d
modelfun_V=@(p)V_fxn(p, t, N0, V0); % vertical
modelfun_both = @(p)mu_V_fxn(p, t, N0);
var_in_mean =  @(p)(1/N).*(V_fxn(p, t, N0, V0)); % vertical



% To find fourth order moment, take parameters, simulate 100 trajectories,
% and find the corresponding v4


var_in_var = @(p)(1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2)));



J= @(phat)((0.5*sum((log(2*pi*var_in_mean(pbxform(phat)))) + (((modelfun_mu(pbxform(phat)) - yfxform(mudata))./sqrt(var_in_mean(pbxform(phat)))).^2))) +...
 (0.5*sum((log(2*pi*var_in_var(pbxform(phat)))) + (((modelfun_V(pbxform(phat)) - yfxform(vardata))./sqrt(var_in_var(pbxform(phat)))).^2))));

objfun_J = @(phat)(J(phat));
J = objfun_J(phat);
Jtrue = real(J);


end

