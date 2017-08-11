% This script applies DEIM CUR with incremental QR to the filtered Incart
% data matrices, tracking which annotations are detected in each file. For
% different CUR stopping tolerances.

% The results are saved as a .mat file named
% [patient_ID_list{i} '_CUR_annotation_tracking.mat']

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create cell array of the different record names
patient_ID_list = cell(75,1);

for id = 1:75
    if id > 9
        patient_ID_list{id} = ['I' num2str(id) 'm'];
    else
        patient_ID_list{id} = ['I0' num2str(id) 'm'];
    end
end

% Apply CUR to all patient files
for i = 1:75
    
    load(['Incart_database/' patient_ID_list{i} '_filtered_data_matrix'])
    
    % Normalize the data matrix corresponding to the Lead II
    data_matrix = data_matrix_beat_normalization(info.data_matrix{2});
    
    
    % CUR tolerances to test; note that these values are divided by 10
    % prior to input into the incremental QR code
    CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
%     CUR_stopping_tol = 5e-4; % Incremental QR tolerance of 5e-5
    
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
    save([patient_ID_list{i} '_CUR_annotation_tracking'],'annotations')
    
    clearvars -except patient_ID_list
    
end