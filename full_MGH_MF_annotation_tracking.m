% This script applies DEIM CUR with incremental QR to the filtered MGH-MF
% data matrices, tracking which annotations are detected in each file. Only
% one CUR stopping tolerance is tested here, but the code is written to
% handle multiple stopping tolerances.

% The results are saved as a .mat file named
% [patient_ID_list{i} '_CUR_annotation_tracking.mat']

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create cell array of the different record names
patient_ID_list = cell(250,1);

for id = 1:250
    if id >= 100
        patient_ID_list{id} = ['mgh' num2str(id) 'm'];
    elseif id > 9 && id < 100
        patient_ID_list{id} = ['mgh0' num2str(id) 'm'];
    else
        patient_ID_list{id} = ['mgh00' num2str(id) 'm'];
    end
end

% Apply CUR to all patient files
for i = 1:250
    
    if i == 61 || i == 127 || i == 230 || i == 235 % These records have no annotations or not Matlab-readable data
        continue
    else
        
        load(['MGH-MF_Waveform_Database/' patient_ID_list{i} '_filtered_data_matrix'])
        
        % Select and normalize the data matrix corresponding to the lead of
        % choice (Lead II is the first choice; otherwise, the first signal
        % without NaNs is analyzed)
        lead2 = 0;
        for j = 1:length(info.signals)
            if strcmpi(info.signals{j},'ECG lead II') || strcmpi(info.signals{j},'ECG leadII') || strcmpi(info.signals{j},'II') || strcmp(info.signals{j},'lead II') || strcmp(info.signals{j},'ECG II') || strcmpi(info.signals{j},'ECG 2') || strcmpi(info.signals{j},'ECG lead 2')
                lead2 = j;
            end
        end
        
        nonan = [];
        test_mat = 1;
        
        while test_mat <= length(info.signals)
            if sum(sum(isnan(info.data_matrix{test_mat}))) > 0
                test_mat = test_mat+1;
            else
                nonan = [nonan test_mat];
                test_mat = test_mat + 1;
            end
        end
        if (~isempty(nonan) && lead2 == 0) || (~isempty(nonan) && ~ismember(lead2,nonan))
            display(['Lead II not available for ' patient_ID_list{i} '. Analyzing first available lead without NaNs, ' info.signals{nonan(1)}])
            data_matrix = data_matrix_beat_normalization(info.data_matrix{nonan(1)});
        elseif ~isempty(nonan) && ismember(lead2,nonan)
            data_matrix = data_matrix_beat_normalization(info.data_matrix{lead2}); % looking at lead II
        else
            display(['Unsuitable data - all ECG signals for ' patient_ID_list{i} 'have NaNs'])
            continue
        end
        
        % CUR tolerances to test; note that these values are divided by 10
        % prior to input into the incremental QR code
%         CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
        CUR_stopping_tol = 5e-4; % Incremental QR tolerance of 5e-5
        
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
        
    end
    clearvars -except patient_ID_list
end