function [ both] = mu_V_fxn( p, tsamp, n0 )
% this function computes the expected value of the variance as a function of time
% for the parameters
% defined parameters
t=tsamp;
b = p(1);
d = p(2);
n_2=(-(b+d)/2*(b-d)).*(n0*exp((b-d).*tsamp)) + ((n0.^2)+((b+d)/(2*(b-d)))*n0).*((exp(2*(b-d).*tsamp)));
mu = n0*exp((b-d)*tsamp);
V =( n0.*(b+d).*t.*exp((b-d)*t));
both = horzcat(mu', V');


end