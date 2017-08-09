% This script loads the different synthetic trial matrices and uses DEIM
% CUR with incremental QR to select a representative subset of beats. The
% precent dimension reduction and number of missed classes are calculated
% and the consolidated results are stored in
% original_DEIM-CUR_synthetic_results_table_with_noise.txt
% with rows corresponding to variability levels and columns corresponding
% to variability types.

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Different variability types
trial_type = cell(5,1);
trial_type{1,1} = 'HR';
trial_type{2,1} = 'shift';
trial_type{3,1} = 'amp';
trial_type{4,1} = 'std';
trial_type{5,1} = 'noise';

% Parameter for CUR with incremental QR
CUR_tol = 1e-4; % this is an threshold of 1e-5 in incremental QR

% Different variability levels
trials = [.01,.02,.05,.1,.2,.3,.5];

% Initialize results vectors
missed_classes_orig = NaN(length(trials),4);
perc_reduction_orig = NaN(length(trials),4);


for j = 1:5 % number of variability types

    for i = 1:length(trials)
        if j == 3 && trials(i) == .5 % a trial matrix was not stored for amplitude variability of 50%
            continue
        end
        
        % Load trial matrices
        filename = [ trial_type{j,1} '_' num2str(trials(i)*100) 'percent'];
        
        load(filename)
        
        % Normalize data (this step may be redundant given that the code
        % currently normalizes the data matrix prior to saving)
        mat = data_matrix_beat_normalization(trial.matrix)'; % transpose the matrix such that each column is a beat
        
        % DEIM CUR with incremental QR
        [C_orig,U_orig,R_orig,p_orig,q_orig] = CURfacQR(mat,CUR_tol);
        
        % Count the class representation
        class_count_orig = synthetic_class_counter(q_orig);
        
        % Determine how many classes were missed
        missed_classes_orig(i,j) = sum(class_count_orig == 0);
        
        % Number of CUR-selected beats
        num_beats_orig = length(q_orig);
        
        % Calculate percent dimension reduction from original number of
        % beats
        perc_reduction_orig(i,j) = (1-(num_beats_orig/size(mat,2)))*100;
       
    end
    
end

% Output results in a table with columns in the order of the variability
% types stored in trial_type
andper = repmat('&',7,1);
endline = repmat('\\',7,1);
table_original_DEIM_results = table(trials'*100, andper,perc_reduction_orig(:,1),andper,missed_classes_orig(:,1),andper,perc_reduction_orig(:,2),andper,missed_classes_orig(:,2),andper,perc_reduction_orig(:,3),andper,missed_classes_orig(:,3),andper,perc_reduction_orig(:,4),andper,missed_classes_orig(:,4),andper,perc_reduction_orig(:,5),andper,missed_classes_orig(:,5),endline);
writetable(table_original_DEIM_results,'original_DEIM-CUR_synthetic_results_table_with_noise','Delimiter',' ','WriteRowNames',true)

        
        