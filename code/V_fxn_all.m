function [ varvec ] = V_fxn_all( p, tin, n0, V0 )
% this function computes the expected value of the variance as a function of time
% for the parameters
% defined parameters

uniqN0=unique(n0);
varvec = [];
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
    
    t = tsamp;


V= N_0(1).*(-(b+d)/(b-d)).*(exp((b-d).*t)) + ((N_0(1).*((b+d)./(b-d)))+V0).*(exp(2.*(b-d).*t));
V = V;
 if time(1)~=0
    V = V(2:end);
 end
 varvec = vertcat(varvec, V);

 end
end