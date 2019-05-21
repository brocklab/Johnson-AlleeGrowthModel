figure; hold off
plot(tlong, Nfakelong,'go')
hold on
plot(tlong(igood), modelfungood(theta),'r*')
plot(tlong(icens), modelfuncens(theta),'b')
%%
int = length(tlong)/length(N0);
for j = 1:length(N0)
    tbigtest(j,:) = tlong(int*(j-1)+1:int*j);
end
%%
ll1 = sum(log(normpdf(yfxform(Nfakelong(igood)),yfxform(modelfungood(params)))))
ll2 =  sum(log(normcdf(yfxform(Nfakelonghigh), yfxform(modelfuncens(params)),sigma)));
ll3 = sum(log(normcdf(yfxform(Nfakelonglow), yfxform(modelfuncens(params)),sigma)))
x = modelfungood(params)

llcomb = (sum(log(normpdf(yfxform(Nfakelong(igood)),yfxform(modelfungood(theta)), sigma)))+...
    sum(log(normcdf(yfxform(Nfakelonghigh), yfxform(modelfuncens(theta)),sigma))));
pbxform(phatbest)

%%
LLG = (sum(log(normpdf(yfxform(Nfakelong(igood)),yfxform(modelfunsinggood(gguess)), sigma)))+...
    sum(log(normcdf(yfxform(Nfakelonghigh), yfxform(modelfunsingcens(gguess)),1))));
%%

num_sims = size(tstate,2);
g = 0.03;
for j =1:num_sims
    N0 = state(1,j);
    tmeas = tstate(:,j);
    N_model(:,j) = singleexpmodel(g, state(1,j),tmeas);
end
vert_length = 10000;
Nmeas = reshape(state, vert_length,1);
N_modellong = reshape(N_model, vert_length,1);
diff = N_modellong - Nmeas;

%% test
diff = fit_single_exp_sim( 0.03, state, tstate);