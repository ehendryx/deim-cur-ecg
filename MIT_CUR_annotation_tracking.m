% This script applies DEIM CUR with incremental QR to the filtered MIT-BIH
% Arrhythmia data matrices, testing different CUR tolerances and tracking
% which annotations are detected in each file.

% This code can be run with and without considering flutter waves. To
% remove flutter waves from consideration, set no_flutter = 1 below.

% The results are saved as a .mat file named
% [patient_ID_list{i} '_CUR_annotation_tracking.mat']
%   with results for file 207m saved as
% '207m_no_flutter_CUR_annotation_tracking'
% when no_flutter is set to 1.

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set this value to 1 to exclude flutter waves
no_flutter = 0;

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

% Apply CUR to all patient files
for i = 1:48
    
    % Load filter data matrices
    load(['MIT_Data/' patient_ID_list{i} '_filtered_data_matrix'])
    
    % Select and normalize the data matrix corresponding to the lead of
    % choice
    if i == 14
        data_matrix = data_matrix_beat_normalization(info.data_matrix2); % patient 114m had V5 before MLII
    else
        data_matrix = data_matrix_beat_normalization(info.data_matrix1); % looking at the MLII lead; for patient 102m and 104m, however, this is V5
    end
    
    % Remove flutter waves from consideration in file 207 if desired
    if no_flutter
        if strcmp(patient_ID_list{i},'207m')
        non_flutter = ~strcmp(info.annotations,'!');
        data_matrix = data_matrix(non_flutter,:);
        info.annotations = info.annotations(non_flutter);
        end
    end
    
    % CUR tolerances to test; note that these values are divided by 10
    % prior to input into the incremental QR code
    CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
    
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
    if no_flutter && strcmp(patient_ID_list{i},'207m')
        save([patient_ID_list{i} '_no_flutter_CUR_annotation_tracking'],'annotations')
    else
        save([patient_ID_list{i} '_CUR_annotation_tracking'],'annotations')
    end
    
    clearvars -except patient_ID_list
    
end