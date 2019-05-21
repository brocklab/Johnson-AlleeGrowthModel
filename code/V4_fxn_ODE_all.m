function  V4vec = V4_fxn_ODE_all( p, n0, V0, tin)

uniqN0=unique(n0);
V4vec = [];
b = p(1);
d = p(2);
for j = 1:length(uniqN0)
    for i = 1:length(tin)
      ind(i)=uniqN0(j)==n0(i);
    end
    time = tin(ind);
    N_0=n0(ind);
            % within each initial cell number group, finds the starting
            % indices
            tstart=find(time<=4);
            numreps =length(tstart);
            for k = 1:numreps
                if k <numreps
                tind = time(tstart(k):tstart(k+1)-1);
                end
                if k == numreps
                    tind = time(tstart(k):end);
                end
            if tind(1)~=0
                tsamp = vertcat(0, tind);
            else
                tsamp = time;
            end
                t = tsamp;
                N0 = N_0(1);
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
        if tind(1)~=0
            V4 = V4(2:end);
        end
        V4vec = vertcat(V4vec, V4);
        end
        
        
       
end
end