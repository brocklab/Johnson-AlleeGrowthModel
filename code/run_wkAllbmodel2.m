function [ Nsamp3D,Nstat, Cstat, musimdatavec, varsimdatavec, timesimdatavec] = run_wkAllbmodel2(paramsAwb, tsamp, Ninit, Nsim)
b = paramsAwb(1);
d = paramsAwb(2);
A = paramsAwb(3);
tau = paramsAwb(4);

musimdatavec = [];
varsimdatavec = [];
timesimdatavec = [];
N0vec = [];
Nsamp = [];
for i = 1:length(Ninit)
%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = Nsim;
num_iters = 2000;
take_offs = 0;
state = zeros(num_iters,num_samps);
tstate = zeros(num_iters,num_samps);
state(1,:) = Ninit(i); % at time 0, number of cells =N
tjump(1, :) = 0; % start at time 0
ct_extinct = 0;
for j = 1:num_samps
    N=Ninit(i);
    N0 = N;
    time(1)= 0;
for k = 2:num_iters
    birth_n  = (b*N-(b-d).*N.*((A+tau)./(N+tau))); % birth 
    if birth_n <0
        birth_n = 0;
    end
    death_n = d*N; % death 
    if N==0
        N=0;
    else
        r = rand;
        if r< (birth_n)/(birth_n+death_n)
        N = N+1;
        end
        if r>= (birth_n)/(birth_n+death_n)
        N = N-1;
        end
    end
    state(k, j) = N;
    % set time step to be proportional
    r2=rand;
    tstep =-log(r2)/(birth_n + death_n);
    if tstep == Inf
        % make tstep same as previous?
        tstep = 1;
    end
    time = time + tstep;
    tstate(k,j) = time;
    
    % If N goes below 0, cells go extinct, N=0 throughout
    if N <= 0
       state(k:end,j) = 0;
    end
   
end
    
    thres = 0;
     if state(end,j) > thres
        take_offs= take_offs +1;
     end 
    ind(j) = state(end,j)>thres;


end


tmin(i) = min(tstate(end, :));
    if tmin(i) < tsamp(end)
        error('Cannot sample, insufficient stochastic time. Increase num_iters')
    end


    for m = 1:num_samps
        tstoch = tstate(:,m);
        Nstoch = state(:,m);

    for n = 1:length(tsamp)
        % find nearest tstate that is less that tsamp
        ind =find(tstoch<=tsamp(n),1,'last');
        tfind = tstoch(ind);
        Nsamp(n,m)=Nstoch(ind);
    end

    end
    mu_data = mean(Nsamp,2);
    n_2_data = mean((Nsamp.^2),2);
    var_data = n_2_data - ((mu_data).^2);
    
    for l=1:length(mu_data)
        if mu_data(l)<0
            mu_data(l)=0;
        end
        if var_data(l)<0
            var_data(l)=0;
        end
    end

    musimdatavec = vertcat(musimdatavec, mu_data);
    varsimdatavec = vertcat(varsimdatavec, var_data);
    timesimdatavec = vertcat(timesimdatavec, tsamp');
    Nsamp3D(:,:,i) = Nsamp;
end
  
% remove t=0 from time vector and musimdatavec and varsimdatavec 
i0 = find(timesimdatavec==0);
timesimdatavec(i0)= [];
musimdatavec(i0)=[];
varsimdatavec(i0) = [];



% Now make Nstat from Nsamp3D
Nstat = [];
for i = 1:length(Ninit)
mu_data = mean(Nsamp3D(:,:,i),2);
n_2_data = mean((Nsamp3D(:,:,i).^2),2);
var_data = n_2_data - ((mu_data).^2);
n_3_data = mean((Nsamp3D(:,:,i).^3),2);
n_4_data = mean((Nsamp3D(:,:,i).^4),2);
var4_data = n_4_data - ((mu_data).^4);
Nstat(:,1,i) = mu_data;
Nstat(:,2,i) = n_2_data;
Nstat(:,3,i) = var_data;
Nstat(:,4,i) = n_3_data;
Nstat(:,5,i) = n_4_data;
Nstat(:,6,i) = var4_data;
end

% Now make Cstat from ODEs

for i = 1:length(Ninit)
    N0 = Ninit(i);

V0 = 0;
C_init(1)=N0;
C_init(2)=N0.^2;
C_init(3)= V0;
C_init(4)=N0.^3;
C_init(5)= N0.^4;
C_init(6) = V0;

f = @(t,C) [((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));  % dN/dt
           2.*C(2).*(b-d) + (b+d).*C(1) - 2.*(C(1).^2).*(b-d).*((A+tau)./(C(1)+tau))-(b-d).*C(1).*((A+tau)./(C(1)+tau));
           2.*C(2).*(b-d) + (b+d).*C(1) - 2.*C(2).*(b-d).*((A+tau)./(C(1)+tau))-(b-d).*C(1).*((A+tau)./(C(1)+tau))+...
           -2.*C(1).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))));
           3.*C(4).*(b-d) - 3.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 3.*C(2).*(b+d) - 3.*C(2).*(b-d).*((A+tau))./(C(1)+tau)+...
           + C(1).*(b-d) - C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn3/dt
           4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) - 6.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
           -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) - C(1).*(b-d).*((A+tau)./(C(1)+tau)); %dn4/dt
            4.*C(5).*(b-d) - 4.*C(5).*(b-d).*((A+tau)./(C(1)+tau)) + 6.*C(4).*(b+d) - 6.*C(4).*(b-d).*((A+tau)./(C(1)+tau)) + 4.*C(2).*(b-d)+...
           -4.*C(2).*(b-d).*((A+tau)./(C(1)+tau)) + C(1).*(b+d) - C(1).*(b-d).*((A+tau)./(C(1)+tau))+...
           -4.*((C(1)).^3).*((b-d)*C(1)*(1-((A+tau)./(C(1)+tau))))]; %d


options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:6);
[t,C]=ode45(f, tsamp,C_init, options);
mu_C= C(:,1);
n2_C= C(:,2);
v2_C = C(:,3);
n3_C = C(:,4);
n4_C = C(:,5);
v4_C=C(:,6);

Cstat(:,1,i) = mu_C;
Cstat(:,2,i) = n2_C;
Cstat(:,3,i) = v2_C;
Cstat(:,4,i) = n3_C;
Cstat(:,5,i) = n4_C;
Cstat(:,6,i) = v4_C;
end
end