function [ dNdt ] = odefunAlleeqdiv(t,N, g,A,carcap, q )
%odefunNs outputs the dNs/dt for a given time, Ns, and Nr, with parameters gs,
%ksr, and krs


dNdt = g.*N.*((N./(A./q))-1).*(1-(N./carcap));

end