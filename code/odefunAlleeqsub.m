function [ dNdt ] = odefunAlleeqsub(t,N, g,A,carcap, q )
%odefunNs outputs the dNs/dt for a given time, Ns, and Nr, with parameters gs,
%ksr, and krs


dNdt = g.*N.*(1-((A-q)./N)).*(1-(N./carcap));

end