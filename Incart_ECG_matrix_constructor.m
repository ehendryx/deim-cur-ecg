% This script reads the *.mat, *.info, and *rr.txt files for each record in
% the Incart Database available on PhysioNet
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

% Reference for plotATM.m:
% Abdala, O., Hislop, H. (2014). "plotATM.m" (Version 1.1) [Computer code].
% Available at https://physionet.org/cgi-bin/atm/ATM (Accessed March 30, 2016.)
% ******************

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

plots = 'no'; % Change this to 'yes' to plot filtered data and peak delineations

% Create cell array of the different record names
patient_ID_list = cell(75,1);

for id = 1:75
    if id > 9
        patient_ID_list{id} = ['I' num2str(id) 'm'];
    else
        patient_ID_list{id} = ['I0' num2str(id) 'm'];
    end
end


% Create data matrix files for each patient ID
for i = 1:75
    
    patient_ID = patient_ID_list{i};
    
    filename = patient_ID;
    load(filename,'-mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % See plotATM.m available on PhysioNet to get sampling interval
    % ("interval"), update "val" in light of signal baseline and gain, as
    % well as signal names to be stored in the cell array "signal"
    
    % Dependency warning
    interval = [];
    val_update = 0; % Change this flag once additional code is added to update val
    signal = [];
    
    if isempty(interval) || val_update == 0 || isempty(signal)
        display('The variables "interval", "val", and "signal" are dependent on code from plotATM.m available on PhysioNet: https://physionet.org/cgi-bin/atm/ATM')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time = 0:interval:(size(val,2)-1)*interval;
    
    % Read R-R interval data
    fileID = fopen([filename 'rr.txt']);
    RRinterval = textscan(fileID,'%f %s %f32 %s %f','HeaderLines',0);
    fclose(fileID);
    
    % Store samples indicating beat start and end
    % ****************************
    % This step should be adjusted based on the format of the *rr.txt file
    beat_start_index = RRinterval{1,1}(2:end);
    
    beat_end_index = RRinterval{1,5}(2:end);
    % ****************************
    
    % Store corresponding annotations
    beat_annotation = cell(1,length(beat_start_index));
    for t = 2:length(RRinterval{1,2})
        beat_annotation{1,t-1} = RRinterval{1,2}{t};
    end
    
    % There are 12 waveforms that correspond to ECGs
    ecg_leads = 1:12;
    
    % Process data from each ECG lead
    for p = 1:length(ecg_leads)
        
        if strcmp(plots,'yes')
            figure
            plot(time,val(p,:))
            hold on
            title(['Raw ECG data for' signal(p)],'fontsize',14,'fontweight','bold')
            
            plot(time(beat_start_index),val(p,beat_start_index),'.k')
            hold on
        end
        
        % Filter data
        filt_ECG = filter_data(time,val(p,:));
        
        % Trim filtered data to avoid edge effects
        trim = round(length(filt_ECG)*.05);
        filt_ECG = filt_ECG(trim:end-trim);
        
        % Isolate beat delineations and annotations for trimmed data
        trim_beat_start_index = beat_start_index(beat_start_index>=trim & beat_end_index<=(length(val(p,:))-trim))-trim+1;
        trim_beat_end_index = beat_end_index(beat_start_index>=trim & beat_end_index<=(length(val(p,:))-trim))-trim+1;
        
        trim_beat_annotation = beat_annotation(beat_start_index>=trim & beat_end_index<=(length(val(p,:))-trim));
        
        
        if strcmp(plots,'yes')
            figure
            plot(time(trim:end-trim),filt_ECG)
            hold on
            plot(time(trim_beat_start_index+trim-1),filt_ECG(trim_beat_start_index),'.k')
            hold on
            plot(time(trim_beat_end_index+trim-1),filt_ECG(trim_beat_end_index),'sm','markersize',8)
            title(['Filtered Data - 5% taken off each end for ' signal(p)],'fontsize',14,'fontweight','bold')
            
            figure
            plot(time,val(p,:))
            hold on
            plot(time(trim_beat_start_index+trim-1),val(p,trim_beat_start_index+trim-1),'.g')
            hold on
            plot(time(trim_beat_end_index+trim-1),val(p,trim_beat_end_index+trim-1),'sm','markersize',8)
        end
        
        if strcmp(plots,'yes')
            figure
        end
        
        % Store data in matrix format
        data_matrix{p} = zeros(length(trim_beat_start_index),150);
        time_matrix = zeros(length(trim_beat_start_index),150);
        for k = 1:length(trim_beat_start_index)
            segment_time = time((trim_beat_start_index(k)+trim-1):(trim_beat_end_index(k)+trim-1));
            segment = filt_ECG(trim_beat_start_index(k):trim_beat_end_index(k));
            
            beat_time = linspace(segment_time(1),segment_time(end),150);
            
            data_matrix{p}(k,:) = interp1(segment_time,segment,beat_time,'linear');
            
            time_matrix(k,:) = beat_time;
            
            if strcmp(plots,'yes')
                plot(time_matrix(k,1),data_matrix{p}(k,1),'*b')
                hold on
                plot(time_matrix(k,end),data_matrix{p}(k,end),'sk')
                hold on
                plot(time_matrix(k,:),data_matrix{p}(k,:),'g')
                hold on
                if k < 10 || ~strcmpi(trim_beat_annotation{k},'N')
                    text(time_matrix(k,1),max(data_matrix{p}(k,:))+.2,trim_beat_annotation{k})
                end
            end
        end
    end
    
    % Store ID, time data, filtered ECG data, annotations, and signal names in info struct
    info.patient_ID = patient_ID;
    info.time_matrix = time_matrix;
    info.data_matrix = data_matrix;
    info.annotations = trim_beat_annotation;
    info.signals = signal(ecg_leads);
    
    % Save info struct
    save([patient_ID '_filtered_data_matrix'],'info')
    
    
%     pause
    close all
    clearvars -except patient_ID_list plots
    
end