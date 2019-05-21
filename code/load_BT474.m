% Load BT-474 all data
% This script is meant to compile a large data set of all BT-474
% trajectories as they come in. This script should generate a large
% structure, grouped by initial cell number, with just the raw initial cell
% number, cell number in time, time, cell number plated, and the date of
% the sort

close all, clear all, clc
% load in each of the data sets and find its length and width
[N1, T1] =xlsread('../data/BT-474_Nseed5_5_30.xls');
sz(1,1:2) = size(N1);
sz(1,2)=sz(1,2)-1;

[N2, T2] =xlsread('../data/BT_474_Nseed2_7_1.xls');
sz(2,1:2) = size(N2);
sz(2,2)=sz(2,2)-1;
sz(:,3)= cumsum(sz(:,2));

[Nr, Tr] =xlsread('../data/BT-474_Nseed_range.xls');
sz(3,1:2) = size(Nr);
sz(3,2)=sz(3,2)-1;
sz(:,3)= cumsum(sz(:,2));
%%
% Make one structure where each entry is a trajectory containing:
% 1. time
% 2. cell  number in time
% 3. initial cell  number
% 4. number seeded
% 5. date of sort
% 6. well of sort
% 7. color

% Write a loop for each data set to give these identifiers

for i = 1:sz(1,2) % size matrix, first row, second column
    BT(i).time = round(N1(:,1),0);
    BT(i).cellnum = N1(:,i+1);
    BT(i).N0 = N1(1, i+1);
    BT(i).Nseed = 5; % this will have to be updated on each loop
    BT(i).date = '5-30-18';
    BT(i).well = T1(1, i+1);
end
% this is where you will add more data sets by using the sz vector to
% concatenate them to the BT structure

%% This chunk adds in data from 7-1-18 data set 
for i = sz(1,3)+1:sz(2,3)
    BT(i).time = round(N2(:,1));
    BT(i).cellnum = N2(:,i-sz(1,3)+1);
    BT(i).N0 = N2(1, i-sz(1,3)+1);
    BT(i).Nseed = 2;
    BT(i).date = '7-1-18';
    BT(i).well = T2(1, i-sz(1,3)+1);

end

%% This chunk adds in wide range data from the pilot study 
%( don't run initially, we already have this separately analyzed)
for i = sz(2,3)+1:sz(3,3) % size matrix, first row, second column
    BT(i).time = round(Nr(:,1));
    BT(i).cellnum = Nr(:,i-sz(2,3)+1);
%     BT(i).N0 = N2(1, i-sz(1,3)+1);
    BT(i).Nseed = []; % this will have to be updated on each loop
    BT(i).date = '4-09-18';
    BT(i).well = Tr(1, i-sz(2,3)+1);
end
%%
for i = sz(2,3)+1:sz(3,3)
    cell4 = { 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C8', 'C9', 'C10', 'C11'};
    cell8 = {'C2', 'C3', 'C4', 'C5', 'C6', 'C7'};
    cell16 = {'D2', 'D3', 'D4', 'D5', 'D6', 'D7'};
    cell32 = {'D8', 'D9', 'D10', 'D11', 'E10', 'E11'};
    cell64 = {'E6', 'E7', 'E8', 'E9'};
    cell128 = {'E2', 'E3', 'E4', 'E5'};
    cell256 = {'F2', 'F3', 'F4', 'F5'};
    cell512 = {'F6', 'F7', 'F8', 'F9'};
    cell1024 = {'F10', 'F11', 'G10', 'G11'};
    cell2048 = {'G6', 'G7', 'G8', 'G9', 'G10'};
    cell4096 = {'G2', 'G3', 'G4', 'G5'};
    if contains(BT(i).well, cell4 )
        BT(i).Nseed = 4;
    end
    if contains(BT(i).well, cell8)
        BT(i).Nseed = 8;
    end
        if contains(BT(i).well, cell16)
        BT(i).Nseed = 16;
        end
        if contains(BT(i).well, cell32)
        BT(i).Nseed = 32;
        end
        if contains(BT(i).well, cell64)
        BT(i).Nseed = 64;
        end
        if contains(BT(i).well, cell128)
        BT(i).Nseed = 128;
        end
        if contains(BT(i).well, cell256)
        BT(i).Nseed = 256;
        end
        if contains(BT(i).well, cell512)
        BT(i).Nseed = 512;
        end
        if contains(BT(i).well, cell1024)
        BT(i).Nseed = 1024;
        end
        if contains(BT(i).well, cell2048)
        BT(i).Nseed = 2048;
        end
        if contains(BT(i).well, cell4096)
        BT(i).Nseed = 4096;
        end
        BT(i).N0 = BT(i).Nseed;
end
%% Order by N0
Afields = fieldnames(BT);
Acell = struct2cell(BT);
szarray= size(Acell);            % Notice that the this is a 3 dimensional array.
                            % For MxN structure array with P fields, the size
                            % of the converted cell array is PxMxN

% Convert to a matrix
Acell = reshape(Acell, szarray(1), []);      % Px(MxN)

% Make each field a column
Acell = Acell';                         % (MxN)xP

% Sort by 3rd field "N0"
Acell = sortrows(Acell, 3);

% Put back into original cell array format
Acell = reshape(Acell', szarray);

% Convert to Struct
BT = cell2struct(Acell, Afields, 1);
%%  Add color by N0
for j = 1:length(BT)
    N0(j) = BT(j).N0;
end
colorsets = varycolor(length(unique(N0)));
uniqN0= unique(N0);

for i = 1:length(BT)
    BT(i).color = [];
    for j = 1:length(uniqN0)
        if BT(i).N0==uniqN0(j)
            BT(i).color =colorsets(j,:);
        end
    end
end
%% Save the raw data structure
% this saves the raw data structure for all imported BT-474 data sets
save('../out/BTall.mat', 'BT')
