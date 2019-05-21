function  V4 = V4_fxn_ODE( p, N0, V0, tin)

b= p(1);
d= p(2);
if tin(1)~=0
    t= horzcat(0,tin);
else
    t=tin;
end
%tsamp = tstart:tint:tstart+100;

 %Expected values for <n>, <n^2>, and variance using numerical solvers

C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3)= N0.^3;
C_init(4)= N0.^4;
C_init(5) = V0;
C_init(6)= V0;
f = @(t,C) [(b-d)*C(1); % dn/dt
            2*C(2)*(b-d) + C(1)*(b+d);   % dn2/dt
            3*C(3)*(b-d) + 3*C(2)*(b+d) + C(1)*(b-d); %dn3/dt
            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d); % dn4/dt
            2*C(5)*(b-d) + (b+d)*C(1); %dV2dt
            4*C(4)*(b-d)+ 6*C(3)*(b+d) + 4*C(2)*(b-d) + C(1)*(b+d)-4*(C(1).^3)*((b-d)*C(1))]; %dV4dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, t,C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
n3_C = C(:,3);
n4_C = C(:,4);
v2_C = C(:,5);
v4_C=C(:,6);

V4 = v4_C;
if tin(1)~=0
    V4 = V4(2:end);
end
end