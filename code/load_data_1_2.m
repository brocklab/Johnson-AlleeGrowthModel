%% load data and create structured array for small cell number data
close all, clear all, clc
[N1, T1] =xlsread('../data/KJ1cellwell1_17.xls');
[N2, T2] = xlsread('../data/KJ1and2cellwell1_17');
[N3, T3] = xlsread('../data/KJ2cellwell1_17.xls');

[Ni, Ti] = xlsread('../data/KJ1to2cellwellinit_nums1_17.xls');
% N 1st row= initial cell number first plate
% N 2nd row = number of colonies at last time point first plate
% 3rd & 4th row => same for 1 & 2 cell plate
% 5th & 6th row => same for 2 cell plate
% 7th & 8th row => same for 8 cell plate

% N is offset from T because first column doesn't have header
offset = 1;
col1 = @(s)(find(strcmp(s,T1(1,:))));
col2 = @(s)(find(strcmp(s,T2(1,:))));
col3 = @(s)(find(strcmp(s,T3(1,:))));
r1 = @(p)(find(strcmp(p,T1(2:end,col1('text')))));
r2 = @(p)(find(strcmp(p,T2(2:end,col2('text')))));
r3 = @(p)(find(strcmp(p,T3(2:end,col3('text')))));

%rows=(find(~isnan(N(:,1))));
rows1 = r1('Y');
rows2 = r2('Y');
rows3 = r3('Y');

N1=N1(rows1,:);
N2=N2(rows2,:);
N3=N3(rows3,:);

T1=T1(1:rows1(end)+1,:);
T2=T2(1:rows2(end)+1,:);
T3=T3(1:rows3(end)+1,:);


l1 = size(N1,2)-1; % number of wells in data set 1
l2 = size(N2,2)-1;% number of wells in data set 2
l3 = size(N3,2)-1;% number of wells in data set 2

%% Load in first data set
colorset = varycolor(12);
time = N1(:,1);
for i = 1:l1
    green(i).cellnum = N1(:,i+1);
    green(i).time = time;
    green(i).N0meas = N1(1, i+1);
    green(i).N0actual = [];
    green(i).well = T1(1, i+2);
    green(i).Nseed = 1;
    green(i).cc = [1 0 0]; % red for 1 cell wells
    green(i).plate ='12_20_18-1';
end
% Repeat for dataset 2
time = N2(:,1);
for i = l1+1:l1+l2
    green(i).cellnum = N2(:,i+1-l1);
    green(i).time = time;
    green(i).N0actual = [];
    green(i).N0meas = N2(1, i+1-l1);
    green(i).well = T2(1, i+2-l1);
    green(i).plate ='12_20_18-2';
end
for i = l1+1:l1+l2
    cell1 = { 'B', 'C', 'D'};
    cell2 = {'E','F','G'};
    if contains(green(i).well, cell1 )
        green(i).Nseed = 1;
        green(i).cc = [1 0 0]; % red
    else
         green(i).Nseed = 2;
         green(i).cc = [ 0 1 1]; % green
    end
end
%Repeat for dataset 3
time = N3(:,1);
for i = l1+ l2+1:l1+l2+l3
    green(i).cellnum = N3(:,i+1-l1-l2);
    green(i).time = time;
    green(i).N0actual = [];
    green(i).N0meas = N3(1, i+1-l1-l2);
    green(i).well = T3(1, i+2-l1-l3);
     green(i).Nseed = 2;
    green(i).cc = [0 1 1]; % green for 2 cell wells
    green(i).plate ='12_20_18-3';
end
%% Add in initial cell number and final colonies formed for the 8 cell wells
for j = 1:l1
    green(j).N0actual = Ni(1,j+1);
    green(j).num_colonies = Ni(2,j+1);
end

for j = l1+1:l1+l1
    green(j).N0actual = Ni(3,j-l1+1);
    green(j).num_colonies = Ni(4,j-l1+1);
end

for j = l1+l2+1:l1+l2+l3
    green(j).N0actual = Ni(5,j-l1-l2+1);
    green(j).num_colonies = Ni(6,j-l1-l2+1);
end
%% Add in 8 cell well from 9-28
S = load('../out/green8.mat');
green8= S.green8;
%% Concatenate green8 to green
green12 = green;
ltot = length(green12);
for j = 1:length(green8)
    green(ltot +j) = green8(j);
end
%% Color by N0actual
for j = 1:length(green)
    for i = 1:9
    if green(j).N0actual == i-1
        green(j).color = colorset(i, :);
    end
    end
%     if green(j).N0actual == 1
%         green(j).color = [ 1 0 0]; % red
%     end
%     if green(j).N0actual == 2
%         green(j).color = [ 0 1 1]; % green
%     end
%     if green(j).N0actual == 3
%         green(j).color = [ 1 1 0]; % purple?
%     end
end

%% Arrange by N0 actual...

Afields = fieldnames(green);
Acell = struct2cell(green);
sz = size(Acell);            % Notice that the this is a 3 dimensional array.
                            % For MxN structure array with P fields, the size
                            % of the converted cell array is PxMxN

% Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)

% Make each field a column
Acell = Acell';                         % (MxN)xP

% Sort by 4th field "Nseed"
Acell = sortrows(Acell, 4); % 4th column is the N0 actual column in this case

% Put back into original cell array format
Acell = reshape(Acell', sz);

% Convert to Struct
green = cell2struct(Acell, Afields, 1);


%%
save('../out/green128.mat', 'green')