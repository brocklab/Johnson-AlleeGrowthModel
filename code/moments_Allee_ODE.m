function [ C ] = moments_Allee_ODE( p, N0, V0, t )
% This function takes a set of Alee parameters b, d, and A and outputs the
% 1-4 moments, and the 2nd and 4th order variance as columns in the matrix
% C sampled at each time t

b= p(1);
d= p(2);
A = p(3);

C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3) = V0;
C_init(4)= N0.^3;
C_init(5)= N0.^4;
C_init(6) = V0;


f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
             2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A);   % dn2/dt
            (2*C(2).*(b-d) - (2*C(1).*(b-d)*A) + (C(1).*(b+d))-((b-d)*A)-(2*C(1).*((b-d)*C(1)-(b-d)*A))); % dV/dt
            (3*C(4).*(b-d)) + (3*C(2).*(b+d)-3*C(2)*(b-d)*A) - (3*C(1).*(b-d)*A)-((b-d)*A);%dn3dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4)*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)-((b-d)*A); %dn4dt
            (4*C(5).*(b-d)) + (6*C(4).*(b+d))-(4*C(4).*(b-d)*A)+(4*C(2).*(b-d))-(6*C(2).*(b-d)*A)+ (C(1).*(b+d))+...
            (4*C(1).*(b-d)*A)-((b-d)*A)-(4.*((C(1)).^3)*(((b-d).*C(1)-(b-d)*A)))];

options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, t, C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
v2_C=C(:,3);
n3_C = C(:,4);
n4_C = C(:,5);
v4_C=C(:,6);

end

