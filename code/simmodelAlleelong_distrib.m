function [Nmodellong] = simmodelAlleelong_distrib(p,sigma, tlong, N0)
g = p(1);
A = p(2);
int = length(tlong)/length(N0);
for j = 1:length(N0)
    tbig(j,:) = tlong(int*(j-1)+1:int*j);
end

for i = 1:length(N0)
    gpick = normrnd(g, sigma); % picks a new growth rate for each initial cell number
Nmodel(:,i)= (N0(i)-A).*exp(gpick*tbig(i,:)) + A;
for j = 1:length(Nmodel(:,i))
    if Nmodel(j,i) <=0
        Nmodel(j,i) = 0;
    end
end
end

Nmodellong = reshape(Nmodel,[size(Nmodel,1)*size(Nmodel,2) ,1]);
end