%% load data and create structured array for small cell number data for 2 and 5 cells for 2 different plates
close all, clear all, clc
[N1, T1] =xlsread('../data/data2&5plate1.xls');
[N2, T2] = xlsread('../data/data2&5KJplate.xls');
%%
% N is offset from T because first column doesn't have header
offset = 1;
col1 = @(s)(find(strcmp(s,T1(1,:))));
col2 = @(s)(find(strcmp(s,T2(1,:))));
r1 = @(p)(find(strcmp(p,T1(2:end,col1('text')))));
r2 = @(p)(find(strcmp(p,T2(2:end,col2('text')))));
%rows=(find(~isnan(N(:,1))));
rows1 = r1('Y');
rows2 = r2('Y');
N1=N1(rows1,:);
N2=N2(rows2,:);
T1=T1(1:rows1(end)+1,:);
T2=T2(1:rows2(end)+1,:);
%%
time = N1(:,1);

for i = 1:size(N1,2)-1
    greenc(i).cellnum = N1(:,i+1);
    greenc(i).time = time;
    greenc(i).N0meas = N1(1, i+1);
    greenc(i).well = T1(1, i+2);
end

for i = 1:2
    greenc(i).Nseed = 2;
end

for i = 3:6 
    greenc(i).Nseed = 5;
end
%%
time = N2(:,1);

for i = 7:6+size(N2,2)-1
    greenc(i).cellnum = N2(:,i-5);
    greenc(i).time = time;
    greenc(i).N0meas = N2(1, i-5);
    greenc(i).well = T2(1, i+2-6);
end
%%
for i = 7:12
    greenc(i).Nseed = 2;
end

for i = 13:16
    greenc(i).Nseed = 5;
end
%%
save('../out/greenc.mat', 'greenc')