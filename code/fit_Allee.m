function [ err ] = fit_Allee( params, Nmeas, Nmeas0, tmeas, num_samps)
%This function returns a cost function to me minimized
g = params(1);
A = params(2);
carcap = params(3);

Nmeas0(end+1) = 0;
Nmeas(end+1) = 0;
% Need run through a loop for each sample to produce curve

for i = 1:num_samps
    % isolate each sample to be run through model
    igood = ismember(Nmeas,Nmeas0(i));
    iend = ismember(Nmeas,Nmeas0(i+1))-1;
    ifit = igood:1:iend;
    t = tmeas(ifit);
    Nmeas = Nmeas(ifit);
    [tout,Nmodel(:,i)] = ode45(@(t,Nmodel) odefunAllee(t, Nmodel,g, A, carcap ),t, Nmeas0(i));
    
end



% now reshape both of them into vector form
  N_modelnew = reshape(Nmodel_tot,[length(Nmeas)-1 1]);

 
 err = N_modelnew - Nmeas;



end