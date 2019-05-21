function [Nmodellong] = simmodelsingexplong(g, tlong, igood, N0)


%int = length(tlong)/length(N0);
ind = find(tlong ==0);


    
Nmodel = [];
for i = 1:length(N0)-1
Nmodel_ind= ((N0(i)).*exp(g*tlong(ind(i):ind(i+1)-1)));
Nmodel = vertcat(Nmodel, Nmodel_ind);
end

Nmodel_ind = ((N0(end)).*exp(g*tlong(ind(end):end)));
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