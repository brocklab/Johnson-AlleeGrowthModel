function [ muvec ] = mu_fxn_all( p, tin, n0 )
% this function computes the expected value of <n> as a function of time
% for the parameters
% defined parameters
% No is a vector of length tin that denotes which initial condition was set
% tin is the total time vector for all different N0

uniqN0=unique(n0);
muvec = [];
b = p(1);
d = p(2);
for j = 1:length(uniqN0)
    for i = 1:length(tin)
      ind(i)=n0(i)==uniqN0(j);
    end
    time = tin(ind);
    N_0=n0(ind);
if time(1)~=0
    tsamp = vertcat(0, time);
else
    tsamp = time;
end
mu = N_0(1).*exp((b-d)*tsamp);

if time(1)~=0
    mu = mu(2:end);
end
muvec = vertcat(muvec,mu);
end
end
