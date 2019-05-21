function [V] = V_fxnA(p, tin, N0, V0)
C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3) = V0;
if tin(1)~=0
    tsamp = horzcat(0,tin);
else
    tsamp= tin;
end
b = p(1);
d = p(2);
A = p(3);

f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
             2*C(2)*(b-d) - (2*C(1)*(b-d)*A) + (C(1)*(b+d))-((b-d)*A);   % dn2/dt
            (2*C(2)*(b-d) - (2*C(1)*(b-d)*A) + (C(1)*(b+d))-((b-d)*A)-(2*C(1)*((b-d)*C(1)-(b-d)*A)))]; % dV/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:3);
[t,C]=ode45(f, tsamp,C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
v2_C=C(:,3);

if tin(1)~=0
    V = v2_C(2:end);
else
    V = v2_C;
end
end