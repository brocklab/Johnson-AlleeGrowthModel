%% load data and create structured array for small cell number data
close all, clear all, clc
[N1, T1] =xlsread('../data/smallcellnumberpure1_4.xls');
[N2, T2] = xlsread('../data/smallcellnumberpure8_64.xls');
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
time1 = N1(:,1);
time2 = N2(:,1);
for i = 1:length(N1)-1
    pure(i).cellnum = N1(:,i+1);
    pure(i).time = time1;
    pure(i).N0 = N1(1:2, i+1);
    pure(i).N0avg = round(mean(pure(i).N0));
end

for i = length(N1):(length(N1)+size(N2,2)-2)
    j = i-length(N1);
    pure(i).cellnum = N2(:,j+2);
    pure(i).time = time2;
    pure(i).N0 = N2(1:2, j+2);
    pure(i).N0avg = round(mean(pure(i).N0));
end
%% Label data by FACS attempt plate number
for j = 1:100
    pure(j).Nplate = 1;
    pure(j).cc = [1 0 0]; % red
end
for j = 101:180
    pure(j).Nplate = 2;
    pure(j). cc = [1 0 1]; % magenta
end
for j = 181: 240
    pure(j).Nplate = 4;
    pure(j).cc = [1 1 0]; % yellow
end

for j = 241:270
    pure(j).Nplate = 8;
    pure(j).cc = [ 0 1 0]; % green
end

for j = 271:300
    pure(j).Nplate = 64;
    pure(j).cc = [ 0 0 1];
end
%% Smooth data usin various methods

for j = 1:length(pure)
    raw = pure(j).cellnum;
    time = pure(j).time;
    post_sg = smooth(raw,30,'sgolay');
    post_rl = smooth(raw, 'rlowess');
    post_rl10 = smooth(raw,10, 'rlowess');
    pure(j).post_sg = round(post_sg,0);
    pure(j).post_rl = post_rl;
    pure(j).post_rl10 = post_rl10;
    t = round(time,0);
    days = 0:48:t(end);% length of time for calculating per capita g
    ind = find(ismember(t, days));
        for k = 1:length(ind)-1
            pure(j).percapitag(k) = (post_sg(ind(k+1))-post_sg(ind(k)))./post_sg(ind(k));
        end
   
    
end
%%
save('../out/pure.mat', 'pure')