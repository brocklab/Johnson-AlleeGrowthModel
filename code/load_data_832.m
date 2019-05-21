%% load data and create structured array for small cell number data
close all, clear all, clc
[N1, T1] =xlsread('../data/KJ928expt832cells.xls');
[N2, T2] = xlsread('../data/KJ121expt1632cells.xls');
[N, T] = xlsread('../data/KJ928expt8cellwellinitialnums.xls');
% N 1st row= initial cell number
% N 2nd row = number of colonies at last time point 

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

l1 = size(N1,2)-1; % number of wells in data set 1
l2 = size(N2,2)-1;% number of wells in data set 2
%% Load in first data set

time = N1(:,1);

for i = 1:l1
    green(i).cellnum = N1(:,i+1);
    green(i).time = time;
    green(i).N0meas = N1(1, i+1);
    green(i).N0actual = [];
    green(i).well = T1(1, i+2);
end

for i = 1:l1
    cell8 = { 'B', 'C', 'D'};
    cell1632 = {'E','F','G'};
    cell32 = {'2', '3', '4', '5', '6'};
    cell16 = {'7', '8', '9', '10', '11'};
    if contains(green(i).well, cell8 )
        green(i).Nseed = 8;
        green(i).cc = [1 0 0]; % red
        green(i).plate = '9-28-1';
    end
    if contains(green(i).well, cell1632)
        if contains(green(i).well, cell16)
            green(i).Nseed = 16;
            green(i).cc = [ 0 1 0]; % green
            green(i).plate = '9-28-1';
        else
            green(i).Nseed = 32;
            green(i).cc = [ 0 0 1]; % blue
            green(i).plate = '9-28-1';
        end
    end
%      t = round(green(i).time,0);
%      hrs = 0:48:t(end);% length of time for calculating per capita g
%      ind = find(ismember(t, hrs));
%         for k = 1:length(ind)-1
%             green(i).percapitag(k) = (green(i).cellnum(ind(k+1))-green(i).cellnum(ind(k)))./green(i).cellnum(ind(k));
%         end
   
end
% Repeat for dataset 2

time = N2(:,1);

for i = l1+1:l1+l2
    green(i).cellnum = N2(:,i+1-l1);
    green(i).time = time;
    green(i).N0meas = N2(1, i+1-l1);
    green(i).N0actual = [];
    green(i).well = T2(1, i+2-l1);
end

for i = l1+1:l1+l2
    cell16 = { 'B', 'C', 'D'};
    cell2432 = {'E','F','G'};
    cell24 = {'2', '3', '4', '5', '6'};
    cell32 = {'7', '8', '9', '10', '11'};
    if contains(green(i).well, cell16 )
        green(i).Nseed = 16;
        green(i).cc = [0 1 0]; % red
        green(i).plate = '9-28-2';
    end
    if contains(green(i).well, cell2432)
        if contains(green(i).well, cell24)
            green(i).Nseed = 24;
            green(i).cc = [ 0 1 1]; % 
            green(i).plate = '9-28-2';
        else
            green(i).Nseed = 32;
            green(i).cc = [ 0 0 1]; % blue
            green(i).plate = '9-28-2';
        end
    end
%      t = round(green(i).time,0);
%      hrs = 0:48:t(end);% length of time for calculating per capita g
%      ind = find(ismember(t, hrs));
%         for k = 1:length(ind)-1
%             green(i).percapitag(k) = (green(i).cellnum(ind(k+1))-green(i).cellnum(ind(k)))./green(i).cellnum(ind(k));
%         end
   
end


%% Order by Nseed
Afields = fieldnames(green);
Acell = struct2cell(green);
sz = size(Acell);            % Notice that the this is a 3 dimensional array.
                            % For MxN structure array with P fields, the size
                            % of the converted cell array is PxMxN

% Convert to a matrix
Acell = reshape(Acell, sz(1), []);      % Px(MxN)

% Make each field a column
Acell = Acell';                         % (MxN)xP

% Sort by 6th field "Nseed"
Acell = sortrows(Acell, 6);

% Put back into original cell array format
Acell = reshape(Acell', sz);

% Convert to Struct
green = cell2struct(Acell, Afields, 1);

%% Add in initial cell number and final colonies formed for the 8 cell wells

num_8s = size(N,2)-1; % should be 30

for j = 1:num_8s
    green(j).N0actual = N(1,j+1);
    green(j).num_colonies = N(2,j+1);
end
%% Make cut off time at t = 165 (shortest data point)
for j = 1:length(green)
     ind=green(j).time <= green(110).time(end) ;
     green(j).time = green(j).time(ind);
     green(j).cellnum = green(j).cellnum(ind);
end
%% Make new structure for exporting
green8= green(1:30);


save('../out/green8.mat', 'green8')
%%
save('../out/green832.mat', 'green')