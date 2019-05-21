function [ err ] = fit_singleexpd( param, Nmeas, N0, tmeas)
%This function returns a cost function to be minimized
    d = param(1);
    N_model = singleexpmodeld(d, N0,tmeas);
    err = N_model - Nmeas;



end