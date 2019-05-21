function [ err ] = fit_singleexp( param, Nmeas, N0, tmeas)
%This function returns a cost function to be minimized
    g = param(1);
    N_model = singleexpmodel(g, N0,tmeas);
    err = N_model - Nmeas;



end