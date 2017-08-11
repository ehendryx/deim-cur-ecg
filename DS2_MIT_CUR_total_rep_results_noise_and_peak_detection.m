% This script reads the *_CUR_annotation_tracking.mat files and summarizes
% the results for the DS2 MIT-BIH Arrhythmia data subset with added noise 
% and different beat delineations. In particular, the percent 
% representation of each annotation is calculated on the full results set 
% and also with the exclusion of annotations that have < 3 representatives 
% in a file.

% The results are saved in the file named
% ['DS2_MIT_CUR_correct_annotation_ident_results_' {noise amount and type} '_noise_' {peak detection method} '_peak_detect']

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

% The CUR tolerances tested in beat selection (recall that these values
% are divided by 10 prior to being input into the incremental QR code)
% CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];
CUR_stopping_tol = 5e-2;

% Loop through beat delineation methods
for pd = 1:2
    % Loop through noise types 
    for ns = 1:4
        % Loop through different noise levels
        for nr = 1:length(SNR)
            
            % If no noise is added, don't loop through noise levels
            if ns == 4 && nr > 1
                continue
            else
                
                % Initialize the set of annotations present in the MIT-BIH Arrhythmia data
                total_annotes{1,1} = 'N';
                add_annote = 2;
                
                % Identify all annotations present in the data set
                patient_beat_count = zeros(1,length(patient_ID_list));
                for i = 1:length(patient_ID_list)
                    
                    % Load tracking results
                    if ns == 4
                        load([ patient_ID_list{i} '_no_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'])
                    else
                        load([ patient_ID_list{i} '_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'])
                    end
                    
                    patient_annotes = annotations.type;
                    
                    % track the total number of beats stored in each filtered data matrix
                    patient_beat_count(i) = sum(annotations.full_distribution);
                    
                    % detect new annotations
                    for t = 1:length(patient_annotes)
                        if ~strcmp(patient_annotes{1,t},total_annotes)
                            total_annotes{add_annote} = patient_annotes{1,t};
                            add_annote = add_annote+1;
                        end
                    end
                    
                    clear patient_annotes annotations
                    
                end
                
                % Initialize cell arrays for representation results
                thresh_annote_rep = cell(1,length(total_annotes));
                thresh_annote_count = cell(1,length(total_annotes));
                thresh_annote_patient_presence = cell(1,length(total_annotes));
                
                annote_rep = cell(1,length(total_annotes));
                annote_count = cell(1,length(total_annotes));
                annote_patient_presence = cell(1,length(total_annotes));
                
                % Determine representation in CUR-selected beats for each annotation
                for j = 1:length(total_annotes)
                    annote_rep{1,j} = NaN(length(CUR_stopping_tol),length(patient_ID_list));
                    annote_count{1,j} = NaN(length(CUR_stopping_tol),length(patient_ID_list));
                    annote_patient_presence{1,j} = zeros(1,length(patient_ID_list));
                    
                    thresh_annote_rep{1,j} = NaN(length(CUR_stopping_tol),length(patient_ID_list));
                    thresh_annote_count{1,j} = NaN(length(CUR_stopping_tol),length(patient_ID_list));
                    thresh_annote_patient_presence{1,j} = zeros(1,length(patient_ID_list));
                    
                    % Look at the annotation's results for each record
                    for i = 1:length(patient_ID_list)
                        
                        % Load tracking results
                        if ns == 4
                            load(['DS2_MIT_noise_peak_CUR_results/' patient_ID_list{i} '_no_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'])
                        else
                            load(['DS2_MIT_noise_peak_CUR_results/' patient_ID_list{i} '_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect_CUR_annotation_tracking'])
                        end
                        
                        patient_annotes = annotations.type;
                        
                        % Determine annotation representation
                        if ~strcmp(annotations.type,total_annotes(j))
                            continue
                        else
                            for k = 1:length(CUR_stopping_tol)
                                annote_rep{1,j}(k,i) = (annotations.CUR_annote_count(k,strcmp(annotations.type,total_annotes(j))) > 0); % logical value of whether or not annotation type is represented in CUR
                                annote_count{1,j}(k,i) = annotations.CUR_annote_count(k,strcmp(annotations.type,total_annotes(j)));
                            end
                            annote_patient_presence{1,j}(i) = annotations.full_distribution(strcmp(annotations.type,total_annotes(j))); % total number of annotation-type present in data
                        end
                        
                        % Determine annotation representation only in cases where 3 or more beats have the annotation
                        if ~strcmp(annotations.type,total_annotes(j))
                            continue
                        else
                            if annotations.full_distribution(strcmp(annotations.type,total_annotes(j))) < 3
                                continue
                            else
                                for k = 1:length(CUR_stopping_tol)
                                    thresh_annote_rep{1,j}(k,i) = (annotations.CUR_annote_count(k,strcmp(annotations.type,total_annotes(j))) > 0); % logical value of whether or not annotation type is represented in CUR
                                    thresh_annote_count{1,j}(k,i) = annotations.CUR_annote_count(k,strcmp(annotations.type,total_annotes(j)));
                                end
                                thresh_annote_patient_presence{1,j}(i) = annotations.full_distribution(strcmp(annotations.type,total_annotes(j))); % total number of annotation-type present in data
                            end
                        end
                    end
                    
                     % Summarize results for the annotation of interest
                    total_annote_patient_presence(1,j) = nansum(annote_patient_presence{j}); % total number of beats with the annotation of interest in the data set
                    total_annote_rep_distrib(:,j) = sum(~isnan(annote_rep{1,j}),2); % how many patients have the annotation
                    total_annote_found(:,j) = nansum(annote_rep{1,j},2); % how many patients have the annotation represented in CUR results
                    total_annote_percent_correct(:,j) = total_annote_found(:,j)./total_annote_rep_distrib(:,j); % percentage of patients with the annotation that also have it represented in CUR results
                    
                    % Same as above, but disregarding patients with less than 3 beats present
                    thresh_total_annote_patient_presence(1,j) = nansum(thresh_annote_patient_presence{j});
                    thresh_annote_rep_distrib(:,j) = sum(~isnan(thresh_annote_rep{1,j}),2);
                    thresh_annote_found(:,j) = nansum(thresh_annote_rep{1,j},2);
                    thresh_annote_percent_correct(:,j) = thresh_annote_found(:,j)./thresh_annote_rep_distrib(:,j);
                      
                end
                
                % Store results summary as a struct
                results.patient_beat_count = patient_beat_count; % number of beats considered for each patient
                results.total_annotations = total_annotes; % the beat annotations present in this data set
                
                results.patient_annote_data_count = annote_patient_presence; % the number of beats with a particular annotation per patient
                results.total_annote_data_count = total_annote_patient_presence; % total number of beats with a particular annotation in the data set
                results.patient_annote_CUR_rep = annote_rep; % whether or not a patient has an annotation represented in CUR results
                results.patient_annote_CUR_count = annote_count; % number of beats with an annotation selected in CUR
                results.total_annote_rep_distrib = total_annote_rep_distrib; % how many patients have an annotation
                results.total_annote_CUR_found = total_annote_found; % how many patients have an annotation represented in CUR results
                results.total_annote_percent_correct = total_annote_percent_correct; % fraction of patients with annotation that have it represented in CUR results
                
                % Same as above, but disregarding patients with less than 3 beats present
                results.thresh_patient_annote_data_count = thresh_annote_patient_presence;
                results.thresh_annote_data_count = thresh_total_annote_patient_presence;
                results.thresh_patient_annote_CUR_rep = thresh_annote_rep;
                results.thresh_patient_annote_CUR_count = thresh_annote_count;
                results.thresh_annote_rep_distrib = thresh_annote_rep_distrib;
                results.thresh_annote_CUR_found = thresh_annote_found;
                results.thresh_annote_percent_correct = thresh_annote_percent_correct;
                
                % Save results
                if ns == 4
                    save(['DS2_MIT_CUR_correct_annotation_ident_results_no_noise_' peak_detect{pd} '_peak_detect'],'results')
                else
                    save(['DS2_MIT_CUR_correct_annotation_ident_results_' SNR_text{nr} 'db_' noise_type{ns} '_noise_' peak_detect{pd} '_peak_detect'],'results')
                end

                clear patient_annotes annotations results patient_beat_count ... 
                    total_annotes annote_patient_presence total_annote_patient_presence ...
                    annote_rep annote_count total_annote_rep_distrib total_annote_found ... 
                    total_annote_percent_correct thresh_annote_patient_presence ...
                    thresh_total_annote_patient_presence thresh_annote_rep thresh_annote_count ...
                    thresh_annote_rep_distrib thresh_annote_found thresh_annote_percent_correct ...
                    add_annote
            end
        end
    end
end
