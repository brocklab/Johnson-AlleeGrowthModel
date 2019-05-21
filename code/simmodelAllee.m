function [Nmodel] = simmodelAllee(p, tbig, N0)
g = p(1);
A = p(2);
for i = 1:length(N0)
Nmodel(:,i)= (N0(i)-A).*exp(g*tbig(i,:)) + A;
for j = 1:length(Nmodel(:,i))
    if Nmodel(j,i) <0
        Nmodel(j,i) = 0;
    end
end
end

end