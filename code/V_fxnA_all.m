function [varvec] = V_fxnA_all(p, tin, n0, V0)
uniqN0=unique(n0);
varvec = [];
b = p(1);
d = p(2);
A = p(3);
for j = 1:length(uniqN0)
    for i = 1:length(tin)
      ind(i)=n0(i)==uniqN0(j);
    end
    time = tin(ind);
    N_0=n0(ind);
            % within each initial cell number group, finds the starting
            % indices
            tstart=find(time<=5);
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



C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3) = V0;

            f = @(t,C) [((b-d)*C(1)-(b-d)*A); % dn/dt
                         2*C(2)*(b-d) - (2*C(1)*(b-d)*A) + (C(1)*(b+d))-((b-d)*A);   % dn2/dt
                        (2*C(2)*(b-d) - (2*C(1)*(b-d)*A) + (C(1)*(b+d))-((b-d)*A)-(2*C(1)*((b-d)*C(1)-(b-d)*A)))]; % dV/dt
            options1 = odeset('Refine',1);  
            options = odeset(options1,'NonNegative',1:3);
            [t,C]=ode45(f, tsamp,C_init, options);
            mu_C= C(:,1);
            n2_C= C(:,2);
            v2_C=C(:,3);

            if tind(1)~=0
             v2_C=v2_C(2:end);
            end
            varvec = vertcat(varvec, v2_C);
            end
            
end            
end