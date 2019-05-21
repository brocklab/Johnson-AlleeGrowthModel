function [Nmodellong] = simmodelAlleelong(p, tlong, igood, N0)
g = p(1);
A = p(2);
int = length(tlong)/length(N0);
% for j = 1:length(N0)
%     tbig(j,:) = tlong(int*(j-1)+1:int*j);
% end

ind = find(tlong ==0);


    
Nmodel = [];
for i = 1:length(N0)-1
Nmodel_ind= ((N0(i)-A).*exp(g*tlong(ind(i):ind(i+1)-1)) + A);
Nmodel = vertcat(Nmodel, Nmodel_ind);
end

Nmodel_ind = ((N0(end)-A).*exp(g*tlong(ind(end):end)) + A);
Nmodel = vertcat(Nmodel, Nmodel_ind);
for j = 1:length(Nmodel)
    if Nmodel(j) <=0
        Nmodel(j) = 0;
    end
end


% Nmodellong = reshape(Nmodel,[size(Nmodel,1)*size(Nmodel,2) ,1]);
Nmodellong = Nmodel;
Nmodellong = Nmodellong(igood);
end