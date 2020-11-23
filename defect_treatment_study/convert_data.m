%%Converts data from the pre-2018 format to the new format

%% %% Set path variables
clear variables; 
importpath='S:\Projects\hBN_defects\MatlabFiles\MatlabDataAnalysisLCB\ImportData';
plotswpspath='S:\Projects\hBN_defects\MatlabFiles\MatlabDataAnalysisLCB\General';
addpath(importpath);
addpath(plotswpspath);

datapath = 'S:\Users\Alex\hbn_scans\defect_creation_study_scans\13\H17\post_sem\Sample 13 H17 polarization'; %path to folder containing data
folder = datapath; %Where to save them, should be the same as the original datapath
HWPpath = fullfile(datapath,['HWP_angles.txt']); %define the data path for the HWP angles.  Should be the same folder as the data....
%% Import Data

basename = '';

%Older Style
WPangles = 0:10:180;
WPangles = WPangles';

%Old Style
WPangles = importdata(HWPpath);  %import the HWP angles

nfiles = length(WPangles);
data_files = cell(1,nfiles);

for ff=1:nfiles
    %Older Style
    data_files{ff} = sprintf('%s%02d',basename,WPangles(ff));

    %Old Style
    data_files{ff} = sprintf('%s%0.1f',basename,WPangles(ff));
end

data_files = cellfun(@(file)fullfile(datapath,[file,'.txt']),data_files,'UniformOutput',0); % add path and extension to each file
dstruc = importAMtextfiles(data_files); % import all data files as sweeps in one data structure
dstruc.name = basename; % Rename using root name (for multiple files)

% Find columns for X Y axis limits 
xmin_sdcol = find(strcmp('X minimum',dstruc.sdheader));
xmax_sdcol = find(strcmp('X maximum',dstruc.sdheader));
ymin_sdcol = find(strcmp('Y minimum',dstruc.sdheader));
ymax_sdcol = find(strcmp('Y maximum',dstruc.sdheader));

xres = (dstruc.sdata(1,xmax_sdcol)-dstruc.sdata(1,xmin_sdcol))/(size(dstruc.data_array{1},1)-1);
yres = (dstruc.sdata(1,ymax_sdcol)-dstruc.sdata(1,ymin_sdcol))/(size(dstruc.data_array{1},2)-1);

% Add SD column for WP angle
dstruc.sdheader = [dstruc.sdheader,'HWP angle (deg)'];
dstruc.sdata = [dstruc.sdata,WPangles];
dstruc.fileinfo(3) = dstruc.fileinfo(3)+1;
pol_sdcol = dstruc.fileinfo(3)+5;

%% Extract all data in 3d array

SwpsToAverage = 1:nfiles;
nswps = length(SwpsToAverage);
npts = dstruc.fileinfo(2);
polvec = dstruc.sdata(SwpsToAverage,pol_sdcol)*pi/180; % polarization vector
all_pldata = zeros(npts,npts,nswps);

for ss=1:nswps
    swp = SwpsToAverage(ss);
    all_pldata(:,:,swp) = fliplr(dstruc.data_array{swp}');
    file = [folder '\' num2str(polvec(ss)*180/pi) '.mat'];
    PL = all_pldata(:, :, swp);
    save(file, 'PL');
end