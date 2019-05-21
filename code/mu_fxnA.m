function [ mu ] = mu_fxnA( p, tin, N0 )
% this function computes the expected value of <n> as a function of time
% for the parameters
% defined parameters
if tin(1)~=0
    tsamp = horzcat(0, tin);
else 
    tsamp = tin;
end
b = p(1);
d = p(2);
A = p(3);
C_init(1)=N0;
f = @(t,C) [((b-d)*C(1)-(b-d)*A)]; % dn/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1);
[t,C]=ode45(f, tsamp,C_init, options);
if tin(1)~=0
    mu= C(2:end,1);
else 
    mu = C;
end
end