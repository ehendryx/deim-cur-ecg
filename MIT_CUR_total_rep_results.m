% This script reads the *_CUR_annotation_tracking.mat files and summarizes
% the results. In particular, the percent representation of each annotation
% is calculated on the full results set and also with the exclusion of
% annotations that have < 3 representatives in a file.

% In this file allows the exclusion of flutter annotations in file 207m if
% desired.

% The results are saved in the file named
% MIT_CUR_correct_annotation_ident_results.mat
%   or
% MIT_CUR_no_flutter_correct_annotation_ident_results
% if flutter waves are excluded from analysis on file 207m

% In addition, tables containing percent representation results are also
% saved in the local directory.

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set value to 1 if flutter beats are to be excluded
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

% The CUR tolerances tested in beat selection (recall that these values
% are divided by 10 prior to being input into the incremental QR code)
CUR_stopping_tol = [.5,.1,5e-2,1e-2,5e-3,1e-3,5e-4,1e-4];

% Initialize the set of annotations present in the MIT-BIH Arrhythmia data
total_annotes{1,1} = 'N';
add_annote = 2;

% Identify all annotations present in the data set
patient_beat_count = zeros(1,length(patient_ID_list));
for i = 1:length(patient_ID_list)
    
    % Load tracking results with or without flutter waves
    if no_flutter && strcmp(patient_ID_list{i},'207m')
        load([patient_ID_list{i} '_no_flutter_CUR_annotation_tracking'])
    else
        load([patient_ID_list{i} '_CUR_annotation_tracking'])
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
    
    clearvars -except patient_ID_list total_annotes add_annote CUR_stopping_tol patient_beat_count
    
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
        
        % Load tracking results with or without flutter waves
        if no_flutter && strcmp(patient_ID_list{i},'207m')
            load([patient_ID_list{i} '_no_flutter_CUR_annotation_tracking'])
        else
            load([patient_ID_list{i} '_CUR_annotation_tracking'])
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
    total_annote_percent_correct(:,j) = (total_annote_found(:,j)./total_annote_rep_distrib(:,j))*100; % percentage of patients with the annotation that also have it represented in CUR results
    
    % Same as above, but disregarding patients with less than 3 beats present
    thresh_total_annote_patient_presence(1,j) = nansum(thresh_annote_patient_presence{j});
    thresh_annote_rep_distrib(:,j) = sum(~isnan(thresh_annote_rep{1,j}),2);
    thresh_annote_found(:,j) = nansum(thresh_annote_rep{1,j},2);
    thresh_annote_percent_correct(:,j) = (thresh_annote_found(:,j)./thresh_annote_rep_distrib(:,j))*100;
    
    
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
results.total_annote_percent_correct = total_annote_percent_correct; % percentage of patients with annotation that have it represented in CUR results

% Same as above, but disregarding patients with less than 3 beats present
results.thresh_patient_annote_data_count = thresh_annote_patient_presence;
results.thresh_annote_data_count = thresh_total_annote_patient_presence;
results.thresh_patient_annote_CUR_rep = thresh_annote_rep;
results.thresh_patient_annote_CUR_count = thresh_annote_count;
results.thresh_annote_rep_distrib = thresh_annote_rep_distrib;
results.thresh_annote_CUR_found = thresh_annote_found;
results.thresh_annote_percent_correct = thresh_annote_percent_correct;

% Save results
if no_flutter
    save('MIT_CUR_no_flutter_correct_annotation_ident_results','results')
else
    save('MIT_CUR_correct_annotation_ident_results','results')
end

% Generate table for CUR tol = 5e-4 (incQR tol = 5e-5)
andper = repmat('&',length(patient_ID_list),1);
backslash = repmat('/',length(patient_ID_list),1);
endline = repmat('\\',length(patient_ID_list),1);
table_result = table(andper,results.patient_annote_CUR_count{1,1}(7,:)',backslash,results.patient_annote_data_count{1,1}',andper,results.patient_annote_CUR_count{1,2}(7,:)',backslash,results.patient_annote_data_count{1,2}',andper,results.patient_annote_CUR_count{1,3}(7,:)',backslash,results.patient_annote_data_count{1,3}',andper,results.patient_annote_CUR_count{1,4}(7,:)',backslash,results.patient_annote_data_count{1,4}',andper,results.patient_annote_CUR_count{1,5}(7,:)',backslash,results.patient_annote_data_count{1,5}',andper,results.patient_annote_CUR_count{1,6}(7,:)',backslash,results.patient_annote_data_count{1,6}',andper,results.patient_annote_CUR_count{1,7}(7,:)',backslash,results.patient_annote_data_count{1,7}',andper,results.patient_annote_CUR_count{1,8}(7,:)',backslash,results.patient_annote_data_count{1,8}',andper,results.patient_annote_CUR_count{1,9}(7,:)',backslash,results.patient_annote_data_count{1,9}',andper,results.patient_annote_CUR_count{1,10}(7,:)',backslash,results.patient_annote_data_count{1,10}',andper,results.patient_annote_CUR_count{1,11}(7,:)',backslash,results.patient_annote_data_count{1,11}',andper,results.patient_annote_CUR_count{1,12}(7,:)',backslash,results.patient_annote_data_count{1,12}',andper,results.patient_annote_CUR_count{1,13}(7,:)',backslash,results.patient_annote_data_count{1,13}',andper,results.patient_annote_CUR_count{1,14}(7,:)',backslash,results.patient_annote_data_count{1,14}',andper,results.patient_annote_CUR_count{1,15}(7,:)',backslash,results.patient_annote_data_count{1,15}',andper,results.patient_annote_CUR_count{1,16}(7,:)',backslash,results.patient_annote_data_count{1,16}',endline,'RowNames',patient_ID_list);
if no_flutter
    writetable(table_result,'patient_CUR_annote_results_incQRtol_5e-5_no_flutter','Delimiter',' ','WriteRowNames',true)
else
    writetable(table_result,'patient_CUR_annote_results_incQRtol_5e-5','Delimiter',' ','WriteRowNames',true)
end

% Generate tables showing results of annotation percent representation in CUR beat selection
andper = repmat('&',8,1);
endline = repmat('\\',8,1);
% Results disregarding patients with less than 3 beats present
table_thresh_results = table(andper,results.thresh_annote_percent_correct(:,1),andper,results.thresh_annote_percent_correct(:,2),andper,results.thresh_annote_percent_correct(:,3),andper,results.thresh_annote_percent_correct(:,4),andper,results.thresh_annote_percent_correct(:,5),andper,results.thresh_annote_percent_correct(:,6),andper,results.thresh_annote_percent_correct(:,7),andper,results.thresh_annote_percent_correct(:,8),andper,results.thresh_annote_percent_correct(:,9),andper,results.thresh_annote_percent_correct(:,10),andper,results.thresh_annote_percent_correct(:,11),andper,results.thresh_annote_percent_correct(:,12),andper,results.thresh_annote_percent_correct(:,13),andper,results.thresh_annote_percent_correct(:,14),andper,results.thresh_annote_percent_correct(:,15),andper,results.thresh_annote_percent_correct(:,16),endline);
if no_flutter
    writetable(table_thresh_results,'CUR_percent_correct_results_with_count_threshold_no_flutter','Delimiter',' ','WriteRowNames',true)
else
    writetable(table_thresh_results,'CUR_percent_correct_results_with_count_threshold','Delimiter',' ','WriteRowNames',true)
end

% Full results
table_full_results = table(andper,results.total_annote_percent_correct(:,1),andper,results.total_annote_percent_correct(:,2),andper,results.total_annote_percent_correct(:,3),andper,results.total_annote_percent_correct(:,4),andper,results.total_annote_percent_correct(:,5),andper,results.total_annote_percent_correct(:,6),andper,results.total_annote_percent_correct(:,7),andper,results.total_annote_percent_correct(:,8),andper,results.total_annote_percent_correct(:,9),andper,results.total_annote_percent_correct(:,10),andper,results.total_annote_percent_correct(:,11),andper,results.total_annote_percent_correct(:,12),andper,results.total_annote_percent_correct(:,13),andper,results.total_annote_percent_correct(:,14),andper,results.total_annote_percent_correct(:,15),andper,results.total_annote_percent_correct(:,16),endline);
if no_flutter
    writetable(table_full_results,'CUR_percent_correct_results_no_flutter','Delimiter',' ','WriteRowNames',true)
else
    writetable(table_full_results,'CUR_percent_correct_results','Delimiter',' ','WriteRowNames',true)
end