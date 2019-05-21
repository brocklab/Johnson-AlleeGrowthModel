% Load BT-474 all data
% This script is meant to compile a large data set of all BT-474
% trajectories as they come in. This script should generate a large
% structure, grouped by initial cell number, with just the raw initial cell
% number, cell number in time, time, cell number plated, and the date of
% the sort

close all, clear all, clc
%%
% load in each of the data sets and find its length and width
[N1, T1] =xlsread('../data/BT_474_Nseed5_5_30wflags.xls');
sz(1,1:2) = size(N1);
sz(1,2)=sz(1,2)-1;

[N2, T2] =xlsread('../data/BT_474_Nseed2_7_1wflags.xls');
sz(2,1:2) = size(N2);
sz(2,2)=sz(2,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N3, T3] =xlsread('../data/BT_474_Nseed10_7_12wflags.xls');
sz(3,1:2) = size(N3);
sz(3,2)=sz(3,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N4, T4] =xlsread('../data/BT_474_Nseed67_30wflags.xls');
sz(4,1:2) = size(N4);
sz(4,2)=sz(4,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N5, T5] =xlsread('../data/BT_474_Nseed15802wflags.xls');
sz(5,1:2) = size(N5);
sz(5,2)=sz(5,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N6, T6] =xlsread('../data/BT-474_Nseed12_8_16wflags.xls');
sz(6,1:2) = size(N6);
sz(6,2)=sz(6,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N7, T7] =xlsread('../data/BT-474_Nseed8_8_16wflags.xls');
sz(7,1:2) = size(N7);
sz(7,2)=sz(7,2)-1;
sz(:,3)= cumsum(sz(:,2));

[N8, T8] =xlsread('../data/BT-474_Nseed_large.xls');
sz(8,1:2) = size(N8);
sz(8,2)=sz(8,2)-1;
sz(:,3)= cumsum(sz(:,2));

% [Nr, Tr] =xlsread('../data/BT-474_Nseed_range.xls');
% sz(3,1:2) = size(Nr);
% sz(3,2)=sz(3,2)-1;
% sz(:,3)= cumsum(sz(:,2));
%%
% Make one structure where each entry is a trajectory containing:
% 1. time
% 2. cell  number in time
% 3. initial cell  number
% 4. number seeded
% 5. date of sort
% 6. well of sort
% 7. color
% 8. bad fit (0 or 1) 1 = badfit, 0 = good fit (flag == 3)
% 9. die off (0 or 1) 1 = die off, 0 = take off (flag == 0)
% 10. persist (0 or 1) 1 = persist, 0 = normal (flag ==2)


% Write a loop for each data set to give these identifiers

for i = 1:sz(1,2) % size matrix, first row, second column
    BT(i).time = round(N1(1:end-1,1),0);
    BT(i).cellnum = N1(1:end-1,i+1);
    BT(i).N0 = N1(1, i+1);
    BT(i).Nseed = 5; % this will have to be updated on each loop
    BT(i).date = '5-30-18';
    BT(i).well = T1(1, i+1);
    if N1(end,i+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    
    if N1(end,i+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N1(end,i+1);
        BT(i).persist = 0;
    end
    if N1(end,i+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N1(end,i+1)==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end   
end
% this is where you will add more data sets by using the sz vector to
% concatenate them to the BT structure

%% This chunk adds in data from 7-1-18 data set 
for i = sz(1,3)+1:sz(2,3)
    BT(i).time = round(N2(1:end-2,1));
    BT(i).cellnum = N2(1:end-2,i-sz(1,3)+1);
    BT(i).N0 = N2(1, i-sz(1,3)+1);
    BT(i).Nseed = 2;
    BT(i).date = '7-1-18';
    BT(i).well = T2(1, i-sz(1,3)+1);
    if N2(end-1,i-sz(1,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N2(end-1,i-sz(1,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N2(end, i-sz(1,3)+1);
        BT(i).persist = 0;
    end
    if N2(end-1,i-sz(1,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N2(end-1,i-sz(1,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end
%% This chunk adds in data from 7-12-18 data set (10 cells per well)
for i = sz(2,3)+1:sz(3,3)
    BT(i).time = round(N3(1:end-2,1));
    BT(i).cellnum = N3(1:end-2,i-sz(2,3)+1);
    BT(i).N0 = N3(1, i-sz(2,3)+1);
    BT(i).Nseed = 10;
    BT(i).date = '7-12-18';
    BT(i).well = T3(1, i-sz(2,3)+1);
    if N3(end,i-sz(2,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N3(end,i-sz(2,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N3(end, i-sz(2,3)+1);
        BT(i).persist = 0;
    end
    if N3(end,i-sz(2,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N3(end,i-sz(2,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end
%% This chunk adds in data from 7-30-18 data set (6 cells per well)
for i = sz(3,3)+1:sz(4,3)
    BT(i).time = round(N4(1:end-1,1));
    BT(i).cellnum = N4(1:end-1,i-sz(3,3)+1);
    BT(i).N0 = N4(1, i-sz(3,3)+1);
    BT(i).Nseed = 6;
    BT(i).date = '7-30-18';
    BT(i).well = T4(1, i-sz(3,3)+1);
    if N4(end,i-sz(3,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N4(end,i-sz(3,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N3(end, i-sz(3,3)+1);
        BT(i).persist = 0;
    end
    if N4(end,i-sz(3,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N4(end,i-sz(3,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end

%% This chunk adds in data from 8-02-18 data set (15 cells per well)
for i = sz(4,3)+1:sz(5,3)
    BT(i).time = round(N5(1:end-1,1));
    BT(i).cellnum = N5(1:end-1,i-sz(4,3)+1);
    BT(i).N0 = N5(1, i-sz(4,3)+1);
    BT(i).Nseed = 15;
    BT(i).date = '8-02-18';
    BT(i).well = T5(1, i-sz(4,3)+1);
    if N5(end,i-sz(4,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N5(end,i-sz(4,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N5(end, i-sz(4,3)+1);
        BT(i).persist = 0;
    end
    if N5(end,i-sz(4,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N5(end,i-sz(4,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end
%% This chunk adds in data from 8-16-18 data set (12 cells per well)
for i = sz(5,3)+1:sz(6,3)
    BT(i).time = round(N6(1:end-1,1));
    BT(i).cellnum = N6(1:end-1,i-sz(5,3)+1);
    BT(i).N0 = N6(1, i-sz(5,3)+1);
    BT(i).Nseed = 12;
    BT(i).date = '8-16-18';
    BT(i).well = T6(1, i-sz(5,3)+1);
    if N6(end,i-sz(5,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N6(end,i-sz(5,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N6(end, i-sz(5,3)+1);
        BT(i).persist = 0;
    end
    if N6(end,i-sz(5,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N6(end,i-sz(5,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end
%% This chunk adds in data from 8-16-18 data set (8 cells per well)
for i = sz(6,3)+1:sz(7,3)
    BT(i).time = round(N7(1:end-1,1));
    BT(i).cellnum = N7(1:end-1,i-sz(6,3)+1);
    BT(i).N0 = N7(1, i-sz(6,3)+1);
    BT(i).Nseed = 8;
    BT(i).date = '8-16-18';
    BT(i).well = T7(1, i-sz(6,3)+1);
    if N7(end,i-sz(6,3)+1) ==3
        BT(i).badfit = 1;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).N0 = NaN;
    end
    
    if N7(end,i-sz(6,3)+1) ==0
        BT(i).badfit = 0;
        BT(i).dieoff = 1;
        BT(i).tdieoff = N7(end, i-sz(6,3)+1);
        BT(i).persist = 0;
    end
    if N7(end,i-sz(6,3)+1) ==2
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 1;
    end
    if N7(end,i-sz(6,3)+1) ==1
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
    end
    

end
%% This chunk adds in data from large data set
for j = 1:length(BT)
    BT(j).V0 = 0;
end

for i = sz(7,3)+1:sz(8,3)
    BT(i).time = round(N8(1:end-1,1));
    BT(i).cellnum = N8(1:end-1,i-sz(7,3)+1);
    BT(i).N0 = N8(1, i-sz(7,3)+1);
    BT(i).date = '8-16-18';
    BT(i).well = T8(1, i-sz(7,3)+1);
        BT(i).badfit = 0;
        BT(i).dieoff = 0;
        BT(i).tdieoff = NaN;
        BT(i).persist = 0;
        BT(i).V0 = 0;

end
for i = sz(7,3)+1:sz(8,3)
    cell512 = { 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'C8', 'C9', 'C10', 'C11',...
        'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7',...
        'D8', 'D9', 'D10', 'D11'};
        cell1024= {'E10', 'E11', 'E6', 'E7', 'E8', 'E9', 'E2', 'E3', 'E4', 'E5','F2', 'F3', 'F4', 'F5'...
       'F6', 'F7', 'F8', 'F9','F10', 'F11', 'G10', 'G11', 'G6', 'G7', 'G8', 'G9', 'G10', 'G2', 'G3', 'G4', 'G5'};
    if contains(BT(i).well, cell512 )
        BT(i).Nseed = 512;
    end
    if contains(BT(i).well, cell1024)
        BT(i).Nseed = 1024;
    end
end

%% DONT RUN This chunk adds in wide range data from the pilot study 
%( don't run initially, we already have this separately analyzed)
for i = sz(2,3)+1:sz(3,3) % size matrix, first row, second column
    BT(i).time = round(Nr(:,1));
    BT(i).cellnum = Nr(:,i-sz(2,3)+1);
%     BT(i).N0 = N2(1, i-sz(1,3)+1);
    BT(i).Nseed = []; % this will have to be updated on each loop
    BT(i).date = '4-09-18';
    BT(i).well = Tr(1, i-sz(2,3)+1);
end
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
