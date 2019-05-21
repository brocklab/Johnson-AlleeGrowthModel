close all, clear all, clc
[N1] =xlsread('../data/popdoublings25.xls');

Nseed = N1(:,1);
numdoublings = N1(:,2);
 dbl2 = [];
 N2 = [];
 
 dbl5 = [];
 N5 = [];
for j = 1:length(Nseed)
    if Nseed(j) == 2
        N2 = vertcat(N2,2);
        dbl2 = vertcat(dbl2, numdoublings(j));
    end
    if Nseed(j) == 5
        N5 = vertcat(N5,2);
        dbl5 = vertcat(dbl5, numdoublings(j));
    end
end
%%
figure;
plot(Nseed, numdoublings, 'o')
xlim([ 1 6])
xlabel ('N_{seed}')
ylabel('Number of Population Doublings')

figure;
boxplot(numdoublings, Nseed)
