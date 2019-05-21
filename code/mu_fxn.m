function [ mu ] = mu_fxn( p, tin, n0 )
% this function computes the expected value of <n> as a function of time
% for the parameters
% defined parameters
b = p(1);
d = p(2);
if tin(1)~=0
    tsamp = horzcat(0, tin);
else
    tsamp = tin;
end
mu = n0*exp((b-d)*tsamp);
mu = mu';

if tin(1)~=0
    mu = mu(2:end);
end


end

