%% load data and create structured array for small cell number data
%specifically for BT-474s
close all, clear all, clc
[N1, T1] =xlsread('../data/BT-474psdataIP.xls');


% N is offset from T because first column doesn't have header
offset = 1;
col1 = @(s)(find(strcmp(s,T1(1,:))));
r1 = @(p)(find(strcmp(p,T1(2:end,col1('text')))));
%rows=(find(~isnan(N(:,1))));
rows1 = r1('Y');
N1=N1(rows1,:);
T1=T1(1:rows1(end)+1,:);

l1 = size(N1,2)-1; % number of wells in data set 1
%% Load in first data set

time = N1(:,1);

for i = 1:l1
    BT(i).cellnum = N1(:,i+1);
    BT(i).time = time;
    BT(i).N0meas = N1(1, i+1);
    BT(i).N0actual = [];
    BT(i).well = T1(1, i+2);
end

for i = 1:l1
    cell4 = { 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C8', 'C9', 'C10', 'C11'};
    cell8 = {'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
    cell16 = {'D2', 'D3', 'D4', 'D5', 'D6', 'D7'};
    cell32 = {'D8', 'D9', 'D10', 'D11', 'E10', 'E11'};
    cell64 = {'E6', 'E7', 'E8', 'E9'};
    cell128 = {'E2', 'E3', 'E4', 'E5'};
    if contains(BT(i).well, cell4 )
        BT(i).Nseed = 4;
        BT(i).cc = [1 0 0]; % red
        BT(i).plate = '4-09-18';
    end
    if contains(BT(i).well, cell8)
        BT(i).Nseed = 8;
        BT(i).cc = [0 1 0]; % green
        BT(i).plate = '4-09-18';
    end
        if contains(BT(i).well, cell16)
        BT(i).Nseed = 16;
        BT(i).cc = [0 0 1]; % green
        BT(i).plate = '4-09-18';
        end
        if contains(BT(i).well, cell32)
        BT(i).Nseed = 32;
        BT(i).cc = [0.5 .5 0]; 
        BT(i).plate = '4-09-18';
        end
        if contains(BT(i).well, cell64)
        BT(i).Nseed = 64;
        BT(i).cc = [0 0.5 0.5]; % green
        BT(i).plate = '4-09-18';
        end
        if contains(BT(i).well, cell128)
        BT(i).Nseed = 128;
        BT(i).cc = [0.5 0 0.5]; % green
        BT(i).plate = '4-09-18';
        end
    
        
    end
%      t = round(green(i).time,0);
%      hrs = 0:48:t(end);% length of time for calculating per capita g
%      ind = find(ismember(t, hrs));
%         for k = 1:length(ind)-1
%             green(i).percapitag(k) = (green(i).cellnum(ind(k+1))-green(i).cellnum(ind(k)))./green(i).cellnum(ind(k));
%         end
   
%% Remove other large cell wells from analysis
for j = 1:length(BT)
    if isempty(BT(j).Nseed)
        keep(j)=false;
    else
        keep(j)=true;
    end
end

BT= BT(keep);

%% Order by Nseed
Afields = fieldnames(BT);
Acell = struct2cell(BT);
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
BT = cell2struct(Acell, Afields, 1);



%% Make new structure for exporting



save('../out/BTps.mat', 'BT')
