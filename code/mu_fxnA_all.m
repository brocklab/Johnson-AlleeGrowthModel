function [ muvec ] = mu_fxnA_all( p, tin, n0 )
% this function computes the expected value of <n> as a function of time
% for the parameters
% defined parameters

uniqN0=unique(n0);
muvec = [];
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
            f = @(t,C) [((b-d)*C(1)-(b-d)*A)]; % dn/dt
            options1 = odeset('Refine',1);  
            options = odeset(options1,'NonNegative',1);
            [t,C]=ode45(f, tsamp,C_init, options);
            mu=C;
            if tind(1)~=0
             mu= mu(2:end);
            end
            muvec = vertcat(muvec, mu);
            end
        
        
       
end
end