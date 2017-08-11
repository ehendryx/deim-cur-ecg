% This script applies DEIM CUR with incremental QR to the DS2 MIT-BIH
% Arrhythmia data subset with different types and levels of added noise and
% different beat delineations. It tests different CUR tolerances, tracking
% which annotations are detected in each file.

% The results are saved as a .mat file named
% [patient_ID_list{i} '_' {noise amount and type} '_noise_' {peak detection method} '_peak_detect_CUR_annotation_tracking']

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% CUR tolerances to test; note that these values are divided by 10
% prior to input into the incremental QR code
% CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
CUR_stopping_tol = 5e-2; % This is multiplied by 1e-10 in the incremental QR algorithm

% Look at annotations results per patient
for i = 1:length(patient_ID_list)
    
    % Loop through different delineation types
    for pd = 1:2
        % Loop through different noise types
        for ns = 1:4
            % Loop through different noise levels
            for nr = 1:length(SNR)
                
                if ns == 4
                    load([ patient_ID_list{i} '_filtered_data_matrix_no_noise_' peak_detect{pd} '_peak_detect'])
                else
                    load([ patient_ID_list{i} '_filtered_data_matrix_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect'])
                end
                
                % Select and normalize the data matrix corresponding to the lead of
                % choice
                if strcmp(patient_ID_list{i},'114m')
                    data_matrix = data_matrix_beat_normalization(info.data_matrix2); % patient 114m had V5 before MLII
                else
                    data_matrix = data_matrix_beat_normalization(info.data_matrix1); % looking at the MLII lead; for patient 102m and 104m, however, this is V5
                end
                
                % Identify the different annotations present in the data
                annotes = unique(info.annotations);
                
                % Determine the number of beats with each annotation
                for j = 1:length(annotes)
                    annote_distrib(j) = sum(strcmp(info.annotations,annotes{1,j}));
                end
                
                % Apply CUR to identify representative beat morphologies
                q = cell(length(CUR_stopping_tol),1);
                for k = 1:length(CUR_stopping_tol)
                    [C,U,R,p,q{1,k}] = CURfacQR(data_matrix',CUR_stopping_tol(k));
                    
                    % Identify corresponding beat annotations
                    selected_q_annotes = cell(1,length(q{1,k}));
                    for t = 1:length(q{1,k})
                        selected_q_annotes{1,t} = info.annotations{q{1,k}(t)};
                    end
                    
                    % Count the number of beats selected for each annotation type
                    for j = 1:length(annotes)
                        annote_count(k,j) = sum(strcmp(selected_q_annotes,annotes{1,j}));
                    end
                end
                
                % Store results in a struct
                annotations.type = annotes; % types of annotations in file
                annotations.full_distribution = annote_distrib; % total number of beats present for each annotation
                annotations.CUR_q = q; % indices of selected beats
                annotations.CUR_annote_count = annote_count; % number of beats selected from each annotation
                
                % Save annotation tracking results
                if ns == 4
                    save([ patient_ID_list{i} '_no_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'],'annotations')
                else
                    save([ patient_ID_list{i} '_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'],'annotations')
                end
                
                clear data_matrix annotes annote_distrib C U R p q selected_q_annotes annote_count annotations
            end
        end
    end
    
end