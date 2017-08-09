% This script reads the *.mat, *.info, and *rr.txt files for each record in
% the MIT-BIH Arrhythmia Database available on PhysioNet
% (https://physionet.org/cgi-bin/atm/ATM), filters each ECG lead in data,
% divides the data into individual beats using the provided annotations,
% and saves the data from each lead in matrix format.
% Note: For simplicity, all of the *.mat, *.info, and *rr.txt files are
% given a uniform naming structure for ease of reading in the data.

% The resulting data matrices, the corresponding time matrix, and the
% corresponding annotations are saved in the file named:
% [patient_ID '_filtered_data_matrix.mat']

% ******************
% **IMPORTANT DEPENDENCY NOTE**:
% Part of this script (denoted below) is dependent on an
% exerpt of the plotATM.m file available on PhysioNet. This code is not
% included here.
% ******************

% Reference for the MIT-BIH Arrhythmia Database:
% Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database.
% IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001). (PMID: 11446209)

% Reference for Physionet:
% Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG,
% Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and
% PhysioNet: Components of a New Research Resource for Complex Physiologic
% Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages;
% http://circ.ahajournals.org/content/101/23/e215.full]; 2000 (June 13).

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

plots = 'yes'; % Change this to 'no' to supress plots

% Create cell array of the different record names
patient_ID_list = cell(48,1);

patient_ID_list{1} = '100m';
patient_ID_list{2} = '101m';
patient_ID_list{3} = '102m';
patient_ID_list{4} = '103m';
patient_ID_list{5} = '104m';
patient_ID_list{6} = '105m';
patient_ID_list{7} = '106m';
patient_ID_list{8} = '107m';
patient_ID_list{9} = '108m';
patient_ID_list{10} = '109m';
patient_ID_list{11} = '111m';
patient_ID_list{12} = '112m';
patient_ID_list{13} = '113m';
patient_ID_list{14} = '114m';
patient_ID_list{15} = '115m';
patient_ID_list{16} = '116m';
patient_ID_list{17} = '117m';
patient_ID_list{18} = '118m';
patient_ID_list{19} = '119m';
patient_ID_list{20} = '121m';
patient_ID_list{21} = '122m';
patient_ID_list{22} = '123m';
patient_ID_list{23} = '124m';
patient_ID_list{24} = '200m';
patient_ID_list{25} = '201m';
patient_ID_list{26} = '202m';
patient_ID_list{27} = '203m';
patient_ID_list{28} = '205m';
patient_ID_list{29} = '207m';
patient_ID_list{30} = '208m';
patient_ID_list{31} = '209m';
patient_ID_list{32} = '210m';
patient_ID_list{33} = '212m';
patient_ID_list{34} = '213m';
patient_ID_list{35} = '214m';
patient_ID_list{36} = '215m';
patient_ID_list{37} = '217m';
patient_ID_list{38} = '219m';
patient_ID_list{39} = '220m';
patient_ID_list{40} = '221m';
patient_ID_list{41} = '222m';
patient_ID_list{42} = '223m';
patient_ID_list{43} = '228m';
patient_ID_list{44} = '230m';
patient_ID_list{45} = '231m';
patient_ID_list{46} = '232m';
patient_ID_list{47} = '233m';
patient_ID_list{48} = '234m';


% Create data matrix files for each patient ID
for pid = 1:48
    
    patient_ID = patient_ID_list{pid};
    
    filename = patient_ID;
    load(filename,'-mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % See plotATM.m available on PhysioNet to get sampling interval
    % ("interval") and to update "val" in light of signal baseline and gain
    
    % Dependency warning
    interval = [];
    val_update = 0; % Change this flag once additional code is added to update val
    
    if isempty(interval) || isempty(val)
        display('The variables "interval" and "val" are dependent on code from plotATM.m available on PhysioNet: https://physionet.org/cgi-bin/atm/ATM')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    time = 0:interval:(size(val,2)-1)*interval;
    
    % Read R-R interval data
    fileID = fopen([filename 'rr.txt']);
    RRinterval = textscan(fileID,'%s %s %f32 %s %s','HeaderLines',0);
    fclose(fileID);
    
    % Convert R-R interval time stamp to sample number 
    
    % ****************************
    % This step should be adjusted based on the format of the *rr.txt file
    beat_start_times = zeros(length(RRinterval{1,1}(2,end)),1); % skip first entry
    for j = 1:length(RRinterval{1,1})-1
        expand = datevec(RRinterval{1,1}(j+1),'MM:SS.FFF'); % output: [year, month, date, hour, min, sec]
        beat_start_times(j) = 3600*expand(4) + 60*expand(5) + expand(6); % convert to seconds
        [close_start,beat_start_index(j)] = min(abs(beat_start_times(j) - time));
    end
    
    beat_end_times = zeros(length(RRinterval{1,5}(2,end)),1);
    for j = 1:length(RRinterval{1,5})-1
        expand = datevec(RRinterval{1,5}(j+1),'MM:SS.FFF'); % output: [year, month, date, hour, min, sec]
        beat_end_times(j) = 3600*expand(4) + 60*expand(5) + expand(6); % convert to seconds
        [close_end,beat_end_index(j)] = min(abs(beat_end_times(j) - time));
    end
    % ****************************
    
    % Store corresponding annotations
    beat_annotation = cell(1,length(beat_start_times));
    for t = 2:length(RRinterval{1,2})
        beat_annotation{1,t-1} = RRinterval{1,2}{t};
    end
    
    if strcmp(plots,'yes')
        figure
        plot(time,val(1,:))
        hold on
        plot(time,val(2,:),'r')
        hold on
        title('Raw ECG data','fontsize',14,'fontweight','bold')
        
        plot(time(beat_start_index),val(1,beat_start_index),'.k')
        hold on
        plot(time(beat_start_index),val(2,beat_start_index),'.g')
        hold on
    end
    
    % Filter data
    filt_ECG1 = filter_data(time,val(1,:));
    filt_ECG2 = filter_data(time,val(2,:));
    
    % Trim filtered data to avoid edge effects
    trim = round(length(filt_ECG1)*.05);
    filt_ECG1 = filt_ECG1(trim:end-trim);
    filt_ECG2 = filt_ECG2(trim:end-trim);
    
    % Isolate beat delineations and annotations for trimmed data
    trim_beat_start_index = beat_start_index(beat_start_index>=trim & beat_end_index<=(length(val(1,:))-trim))-trim+1;
    trim_beat_end_index = beat_end_index(beat_start_index>=trim & beat_end_index<=(length(val(1,:))-trim))-trim+1;
    
    trim_beat_annotation = beat_annotation(beat_start_index>=trim & beat_end_index<=(length(val(1,:))-trim));
    
    
    if strcmp(plots,'yes')
        figure
        plot(time(trim:end-trim),filt_ECG1)
        hold on
        plot(time(trim:end-trim),filt_ECG2,'r')
        hold on
        plot(time(trim_beat_start_index+trim-1),filt_ECG1(trim_beat_start_index),'.k')
        hold on
        plot(time(trim_beat_start_index+trim-1),filt_ECG2(trim_beat_start_index),'.g')
        hold on
        plot(time(trim_beat_end_index+trim-1),filt_ECG1(trim_beat_end_index),'sm','markersize',8)
        title('Filtered Data - 5% taken off each end','fontsize',14,'fontweight','bold')
        
        figure
        plot(time,val(1,:))
        hold on
        plot(time,val(2,:),'r')
        hold on
        plot(time(trim_beat_start_index+trim-1),val(1,trim_beat_start_index+trim-1),'.g')
        hold on
        plot(time(trim_beat_start_index+trim-1),val(2,trim_beat_start_index+trim-1),'.k')
        hold on
        plot(time(trim_beat_end_index+trim-1),val(1,trim_beat_end_index+trim-1),'sm','markersize',8)
    end
    
    if strcmp(plots,'yes')
        figure
    end
    
    % Store data in matrix format
    data_matrix = zeros(length(trim_beat_start_index),150);
    time_matrix = zeros(length(trim_beat_start_index),150);
    for k = 1:length(trim_beat_start_index)
        segment_time = time((trim_beat_start_index(k)+trim-1):(trim_beat_end_index(k)+trim-1));
        segment1 = filt_ECG1(trim_beat_start_index(k):trim_beat_end_index(k));
        segment2 = filt_ECG2(trim_beat_start_index(k):trim_beat_end_index(k));
        
        beat_time = linspace(segment_time(1),segment_time(end),150);
        
        data_matrix1(k,:) = interp1(segment_time,segment1,beat_time,'linear');
        
        data_matrix2(k,:) = interp1(segment_time,segment2,beat_time,'linear');
        
        time_matrix(k,:) = beat_time;
        
        if strcmp(plots,'yes')
            plot(time_matrix(k,1),data_matrix1(k,1),'*b')
            hold on
            plot(time_matrix(k,end),data_matrix1(k,end),'sk')
            hold on
            plot(time_matrix(k,:),data_matrix1(k,:),'g')
            hold on
            
            plot(time_matrix(k,1),data_matrix2(k,1),'*m')
            hold on
            plot(time_matrix(k,end),data_matrix2(k,end),'sc')
            hold on
            plot(time_matrix(k,:),data_matrix2(k,:),'r')
            hold on
        end
    end
    
    % Store ID, time data, filtered ECG data, and annotations in info struct
    info.patient_ID = patient_ID;
    info.time_matrix = time_matrix;
    info.data_matrix1 = data_matrix1;
    info.data_matrix2 = data_matrix2;
    info.annotations = trim_beat_annotation;
    
    % Save info struct
    save([patient_ID '_filtered_data_matrix'],'info')
    
    pause
    close all
    clearvars -except patient_ID_list
    
end