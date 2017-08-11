% This script adds noise the em, ma, and bw noise waveforms from the
% MIT-BIH Noise Stress Test Database to the DS2 subset of the MIT-BIH
% Arrhythmia Database and divides the data in one of two ways: the
% PhysioNet-provided annotations or an automated peak finder.
% The code reads the *.mat, *.info, and *rr.txt files for each record in
% the MIT-BIH Arrhythmia Database available on PhysioNet
% (https://physionet.org/cgi-bin/atm/ATM), filters each ECG lead in data,
% divides the data into individual beats,
% and saves the data from each lead in matrix format.
% Note: For simplicity, all of the *.mat, *.info, and *rr.txt files are
% given a uniform naming structure for ease of reading in the data.

% The resulting data matrices, the corresponding time matrix, and the
% corresponding annotations are saved in the file named:
% [patient_ID '_filtered_data_matrix_' {noise amount and type} '_noise_' {peak detection method} '_peak_detect.mat']

% ******************
% **IMPORTANT DEPENDENCY NOTE**:
% Two portions of this script (denoted below) are dependent on an
% exerpt of the plotATM.m file available on PhysioNet. This code is not
% included here.

% Reference for plotATM.m:
% Abdala, O., Hislop, H. (2014). "plotATM.m" (Version 1.1) [Computer code].
% Available at https://physionet.org/cgi-bin/atm/ATM (Accessed March 30, 2016.)
% ******************

% Reference for the MIT-BIH Noise Stress Test Database:
% G. B. Moody, W. Muldrow, R. G. Mark, A noise stress test for arrhythmia
% detectors, Computers in cardiology 11 (3) (1984) 381?384.

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

plots = 'yes'; % Change this to 'yes' to plot filtered data and peak delineations

% Peak detection options
peak_detect{1} = 'annote';
peak_detect{2} = 'simple';

% Noise levels to test
SNR = [-6 0 6 12 18 24];
SNR_text{1} = 'n6';
SNR_text{2} = '0';
SNR_text{3} = '6';
SNR_text{4} = '12';
SNR_text{5} = '18';
SNR_text{6} = '24';

% Noise types
noise_type{1} = 'em';
noise_type{2} = 'ma';
noise_type{3} = 'bw';
noise_type{4} = 'none';

% Load noise signals from MIT-BIH Noise Stress Test Database
for ns = 1:3
    
    filename = [noise_type{ns} 'm'];
    load(filename,'-mat');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % See plotATM.m available on PhysioNet to get sampling interval
    % ("interval") and to update "val" in light of signal baseline and gain
    
    % Dependency warning
    interval = [];
    val_update = 0; % Change this flag once additional code is added to update val
    
    if isempty(interval) || val_update == 0
        display('The variables "interval" and "val" are dependent on code from plotATM.m available on PhysioNet: https://physionet.org/cgi-bin/atm/ATM')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    noise_signal{ns} = val;
end
noise_signal{4} = zeros(size(val));

% DS2 MIT-BIH set
patient_ID_list = cell(22,1);

patient_ID_list{1} = '100m';
patient_ID_list{2} = '103m';
patient_ID_list{3} = '105m';
patient_ID_list{4} = '111m';
patient_ID_list{5} = '113m';
patient_ID_list{6} = '117m';
patient_ID_list{7} = '121m';
patient_ID_list{8} = '123m';
patient_ID_list{9} = '200m';
patient_ID_list{10} = '202m';
patient_ID_list{11} = '210m';
patient_ID_list{12} = '212m';
patient_ID_list{13} = '213m';
patient_ID_list{14} = '214m';
patient_ID_list{15} = '219m';
patient_ID_list{16} = '221m';
patient_ID_list{17} = '222m';
patient_ID_list{18} = '228m';
patient_ID_list{19} = '231m';
patient_ID_list{20} = '232m';
patient_ID_list{21} = '233m';
patient_ID_list{22} = '234m';

for pid = 1:length(patient_ID_list)
    
    patient_ID = patient_ID_list{pid};
    
    filename = ['../../MIT_Data/' patient_ID];
    load([filename '.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % See plotATM.m available on PhysioNet to get the sampling frequency and sampling interval
    % ("interval") held in the "freqint" array and to update "val" in light of signal baseline and gain
    
    % Dependency warning
    freqint = [];
    interval = [];
    val_update = 0; % Change this flag once additional code is added to update val
    
    if isempty(freqint) || isempty(interval) || val_update == 0
        display('The variables "feqint", "interval" and "val" are dependent on code from plotATM.m available on PhysioNet: https://physionet.org/cgi-bin/atm/ATM')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fs = freqint(1);
    
    time = 0:interval:(size(val,2)-1)*interval;
    
    if strcmp(plots,'yes')
        figure
        plot(time,val(1,:))
        hold on
        plot(time,val(2,:),'r')
        hold on
        title('Raw ECG data','fontsize',14,'fontweight','bold')
    end
    
    %%% Add noise (or no noise) %%%
    for ns = 1:4
        
        % Add this noise for different signal-to-noise (SNR) ratios
        for nr = 1:length(SNR)
            
            % If no noise is added, only one data matrix needs to be
            % constructed per beat delieation method
            if ns == 4 && nr > 1
                continue
            else
                
                % Pad noise, mirroring values at end (symmetrically) to be the same length as clean data
                padded_noise1 = padarray(noise_signal{ns}(1,:),[0 length(val(1,:))-length(noise_signal{ns}(1,:))],'symmetric','post');
                padded_noise2 = padarray(noise_signal{ns}(2,:),[0 length(val(2,:))-length(noise_signal{ns}(2,:))],'symmetric','post');
                
                % Compute power of clean and noisy signals
                Psig1 = rms(val(1,:))^2;
                Psig2 = rms(val(2,:))^2;
                P_noise1 = rms(padded_noise1)^2;
                P_noise2 = rms(padded_noise2)^2;
                
                if ns == 4
                    a1 = 0;
                    a2 = 0;
                else
                    a1 = sqrt(exp(-log(10)*SNR(nr)/10)*Psig1/P_noise1);
                    a2 = sqrt(exp(-log(10)*SNR(nr)/10)*Psig2/P_noise2);
                end
                
                % Add noise to signal with desired SNR
                noisy_sig1 = val(1,:) + a1*padded_noise1;
                noisy_sig2 = val(2,:) + a2*padded_noise2;
                
                % Filter data
                filt_ECG1 = filter_data(time,noisy_sig1);
                filt_ECG2 = filter_data(time,noisy_sig2);
                
                % Trim filtered data to avoid edge effects
                trim = round(length(filt_ECG1)*.05);
                filt_ECG1 = filt_ECG1(trim:end-trim);
                filt_ECG2 = filt_ECG2(trim:end-trim);
                
                % Identify original beat deliations
                fileID = fopen([filename 'rr.txt']);
                RRinterval = textscan(fileID,'%s %s %f32 %s %s','HeaderLines',0);
                fclose(fileID);
                
                % Convert R-R interval time stamp to sample number
                
                % ****************************
                % This step should be adjusted based on the format of the *rr.txt file
                beat_start_times = zeros(length(RRinterval{1,1}(2,end)),1); % skip first entry
                labeled_beat_start_index = zeros(length(RRinterval{1,1}(2,end)),1);
                for j = 1:length(RRinterval{1,1})-1
                    expand = datevec(RRinterval{1,1}(j+1),'MM:SS.FFF'); % output: [year, month, date, hour, min, sec]
                    beat_start_times(j) = 3600*expand(4) + 60*expand(5) + expand(6); % convert to seconds
                    [close_start,labeled_beat_start_index(j)] = min(abs(beat_start_times(j) - time));
                end
                
                beat_end_times = zeros(length(RRinterval{1,5}(2,end)),1);
                labeled_beat_end_index = zeros(length(RRinterval{1,1}(2,end)),1);
                for j = 1:length(RRinterval{1,5})-1
                    expand = datevec(RRinterval{1,5}(j+1),'MM:SS.FFF'); % output: [year, month, date, hour, min, sec]
                    beat_end_times(j) = 3600*expand(4) + 60*expand(5) + expand(6); % convert to seconds
                    [close_end,labeled_beat_end_index(j)] = min(abs(beat_end_times(j) - time));
                end
                % ****************************
                
                % Store corresponding annotations
                full_beat_annotation = cell(1,length(labeled_beat_start_index));
                for t = 2:length(RRinterval{1,2})
                    full_beat_annotation{1,t-1} = RRinterval{1,2}{t};
                end
                
                % Isolate beat delineations for trimmed data
                labeled_trim_beat_start_index = labeled_beat_start_index(labeled_beat_start_index>=trim & labeled_beat_end_index<=(length(val(1,:))-trim))-trim+1;
                labeled_trim_beat_end_index = labeled_beat_end_index(labeled_beat_start_index>=trim & labeled_beat_end_index<=(length(val(1,:))-trim))-trim+1;
                
                % Delineate beats
                for pd = 1:2
                    if strcmp(peak_detect(pd), 'annote')
                        trim_beat_start_index = labeled_trim_beat_start_index;
                        trim_beat_end_index = labeled_trim_beat_end_index;
                        
                        % Isolate beat annotations for trimmed data
                        trim_beat_annotation = full_beat_annotation(labeled_beat_start_index>=trim & labeled_beat_end_index<=(length(val(1,:))-trim));
                        
                        
                    elseif strcmp(peak_detect(pd),'simple')
                        
                        % Find peaks on filtered data
                        if pid == 14
                            peak_indices = find_peaks(time(:),abs(filt_ECG2(:))); % patient 114m had V5 before MLII
                        else
                            peak_indices = find_peaks(time(:),abs(filt_ECG1(:))); % looking at the MLII lead; for patient 102m and 104m, however, this is V5
                        end
                        peak_indices = peak_indices{1};
                        trim_beat_start_index = peak_indices(1:end-1);
                        trim_beat_end_index = peak_indices(2:end);
                        
                        % Assign labels to beats based on nearest annotation
                        % from original data
                        trim_beat_annotation = cell(1,length(trim_beat_start_index));
                        for t = 1:length(trim_beat_start_index)
                            [close_label,close_label_index] = min(abs((trim_beat_start_index(t)+trim-1) - labeled_beat_start_index));
                            trim_beat_annotation{1,t} = RRinterval{1,2}{close_label_index+1};
                        end
                        
                    else
                        display('Invalid Peak Detection Algorithm')
                        break
                    end
                    
                    
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
                        plot(time,noisy_sig1)
                        hold on
                        plot(time,noisy_sig2,'r')
                        hold on
                        plot(time(trim_beat_start_index+trim-1),noisy_sig1(trim_beat_start_index+trim-1),'.g')
                        hold on
                        plot(time(trim_beat_start_index+trim-1),noisy_sig2(trim_beat_start_index+trim-1),'.k')
                        hold on
                        plot(time(trim_beat_end_index+trim-1),noisy_sig1(trim_beat_end_index+trim-1),'sm','markersize',8)
                    end
                    
                    if strcmp(plots,'yes')
                        figure
                    end
                    
                    % Store data in matrix format
                    data_matrix1 = zeros(length(trim_beat_start_index),150);
                    data_matrix2 = zeros(length(trim_beat_start_index),150);
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
                    if ns == 4
                        save([patient_ID '_filtered_data_matrix_no_noise_' peak_detect{pd} '_peak_detect'],'info')
                    else
                        save([patient_ID '_filtered_data_matrix_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect'],'info')
                    end
                    
%                     pause
                    close all
                    clear trim_beat_start_index trim_beat_end_index trim_beat_annotation segment_time segment1 segment2 beat_time time_matrix data_matrix1 data_matrix2 info
                end
                clear padded_noise1 padded_noise2 Psig1 Psig2 P_noise1 P_noise2 a1 a2 noisy_sig1 noisy_sig2 filt_ECG1 filt_ECG2 trim RRinterval ...
                    beat_start_times labeled_beat_start_index beat_end_times labeled_beat_end_index full_beat_annotation labeled_trim_beat_start_index ...
                    labeled_trim_beat_end_index
            end
        end
    end
end