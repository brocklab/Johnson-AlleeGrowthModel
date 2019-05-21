% Profile Likelihood Generator
% Marisa Eisenberg (marisae@umich.edu) - 7/31/16 - updated 6-22-17

function profile = ProfLike(params,profindex,factor, t, data,N0)
% Definitions
%   params = point in parameter space from which to profile (the parameter estimates)
%   profparam = index of the parameter to be profiled
%   factor = the fractional/percent range to profile the parameter over
%   costfun = this is a cost function of only the parameters (the full vector of them)
%     Everythng else (data, ICs, etc.) is fixed for the entire profile so is set outside when 
%     costfun is defined. In ProfLike, we'll put a wrapper on costfun that
%     will fix the profiled parameter.

% Setup
% factor = 0.5;
numpoints = 2;
mudata = data(:,1);
vardata= data(:,2);
N= length(t);
V0 = 0;

pfxform = @(pval)[1 1].*log(pval); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1].*exp(phat); %'backward' parameter transform into model space
modelfun_mu = @(p)mu_fxn(p, t, N0); % <n> for params (p) = b & d, vertical
modelfun_V=@(p)V_fxn(p, t, N0, V0); % vertical
var_in_mean =  @(p)(1/N).*(V_fxn(p, t, N0, V0)); % vertical
var_in_var = @(p)(1/N).*(V4_fxn_ODE(p,N0, V0, t)-(((N-3)./(N-1)).*(V_fxn(p, t, N0, V0).^2)));

J = @(phat) (sum(((mudata-modelfun_mu(pbxform(phat))).^2)./(2.*sqrt(var_in_mean(pbxform(phat)))) + log(sqrt(var_in_mean(pbxform(phat)))) + 0.5*log(2*pi))+...
    sum(((vardata-modelfun_V(pbxform(phat))).^2)./(2.*sqrt(var_in_var(pbxform(phat)))) + log(sqrt(var_in_var(pbxform(phat)))) + 0.5*log(2*pi)));
objfun_J = @(phat)(J(phat));

% Costfun wrapper
if profindex == 1
fixcostfun = @(shortparams,profparamval)...
    objfun_J(pfxform([profparamval,shortparams(1)]));
end

if profindex ==2
    fixcostfun = @(shortparams,profparamval)...
    objfun_J(pfxform([shortparams(1),profparamval(1)]));
end
% Profile
profrangeDown = linspace(params(profindex), params(profindex)*(1-factor),numpoints)'; 
profrangeUp = linspace(params(profindex), params(profindex)*(1+factor),numpoints)';
% split into up and down so we can use last fitted value as starting value for next run
profrange = [profrangeDown profrangeUp];
currfvals = [];
currparams = [];
currflags = [];

for i=1:2
    if profindex ==1 
    paramstemp = params(2);
    end
    if profindex ==2
        paramstemp = params(1);
    end
    for j = 1:numpoints
        [i j]; %track progress
        [phattemp, fvaltemp, flagtemp] = fminsearch(@(p) fixcostfun(p,profrange(j,i)),paramstemp,optimset('MaxFunEvals',5000,'MaxIter',5000));
        paramstemp = pbxform(phattemp);
        currfvals = [currfvals; fvaltemp];
        currflags = [currflags; flagtemp];
        currparams = [currparams; [paramstemp(1:profindex-1),profrange(j,i),paramstemp(profindex:end)]]; %storing the profiled value too, so the output parameter values are easy to run the model with
    end
end

profile = [flipud([profrangeDown currfvals(1:numpoints) currflags(1:numpoints) currparams(1:numpoints,:)]);...
    [profrangeUp currfvals(numpoints+1:end) currflags(numpoints+1:end) currparams(numpoints+1:end,:)]];
end