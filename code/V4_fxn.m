function [V4] = V4_fxn( p, N0, tstart, tint )
% This function simulates 1000 trajectories for the corresponding
% parameters of b and d, and outputs the fourth order moment corresponding
% to that SSA

Ninit = N0;
b= p(1);
d= p(2);


%birth_n = b*N; % birth hazard function
%death_n = d*N; % death hazard function
num_samps = 500;
num_iters = 250;
take_offs = 0;
state = zeros(num_iters,num_samps);
tstate = zeros(num_iters,num_samps);
state(1,:) = Ninit; % at time 0, number of cells =N
tjump(1, :) = 0; % start at time 0
ct_extinct = 0;
for j = 1:num_samps

    N=Ninit;
    N0 = N;
    time(1)= 0;
for k = 2:num_iters
    birth_n = b*N; % birth 
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


% find minimum time "measured"
tmin = min(tstate(end, :));

% Add in t = 0 to your time matrix
tsamp = tstart:tint:tstart+100;
for j = 1:num_samps
    tstoch = tstate(:,j);
    Nstoch = state(:,j);

for i = 1:length(tsamp)
    % find nearest tstate that is less that tsamp
    ind =find(tstoch<=tsamp(i),1,'last');
    tfind = tstoch(ind);
    Nsamp(i,j)=Nstoch(ind);
end

end

mu_data = mean(Nsamp,2);
sum_quad_diff=zeros(length(tsamp),1);
for j = 1:num_samps
    quad_diff = (Nsamp(:,j)-mu_data).^4;
    sum_quad_diff = sum_quad_diff + quad_diff;
end

var_data4 = sum_quad_diff/num_samps;
V4 = var_data4;
mu = mu_data;

end
