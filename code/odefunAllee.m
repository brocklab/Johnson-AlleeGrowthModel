function [ dNdt, percapdNdt ] = odefunAllee(t,N, g,A,carcap )
%odefunNs outputs the dNs/dt for a given time, Ns, and Nr, with parameters gs,
%ksr, and krs


%dNdt = g*N*(1-(A/N))*(1-(N./carcap));
dNdt = g*N*((N/A)-1)*(1-(N./carcap));
percapdDNdt = dNdt/N;

end