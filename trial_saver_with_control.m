% This script calls the function trial_matrix_generatior to generate
% synthetic ECG data with different levels of variability. These matrices are
% saved in the current directory.

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct Control Matrix

% Maximum amount of change allowed
HR_percent_change = 0; % (percent/100)
shift_percent_change = 0;
amp_percent_change = 0;
std_percent_change = 0;
noise_percentage = 0;

control_matrix = trial_matrix_generator(HR_percent_change,shift_percent_change,amp_percent_change,std_percent_change,noise_percentage);

control.matrix = control_matrix;

save('Control','control')


%% Construct Trial Matrices

% The different types of variability
trial_type = cell(5,1);
trial_type{1,1} = 'HR';
trial_type{2,1} = 'shift';
trial_type{3,1} = 'amp';
trial_type{4,1} = 'std';
trial_type{5,1} = 'noise';

% The different levels of variability
trials = [.01,.02,.05,.1,.2,.3,.5];

% Generate a data matrix with each type and level of variability
for j = 1:5 % number of trial_types
    HR_percent_change = 0;
    shift_percent_change = 0;
    amp_percent_change = 0;
    std_percent_change = 0;
    noise_percentage = 0;
    
    for i = 1:length(trials)
        filename = [trial_type{j,1} '_' num2str(trials(i)*100) 'percent'];
        if j == 1
            HR_percent_change = trials(i);
        elseif j == 2
            shift_percent_change = trials(i);
        elseif j == 3
            amp_percent_change = trials(i);
        elseif j == 4
            std_percent_change = trials(i);
        elseif j == 5
            noise_percentage = trials(i);
        end
        
        if j == 3 && trials(i) == .5
            continue
        end
        
        data_matrix = trial_matrix_generator(HR_percent_change,shift_percent_change,amp_percent_change,std_percent_change,noise_percentage);
        
        pause
        close all
        
        % Create a struct to store resulting data matrix
        trial.type = trial_type(j); % type of variability
        trial.percent = trials(i); % level of variability
        trial.matrix = data_matrix; % resulting matrix
        
        save(filename,'trial')
    end
end

