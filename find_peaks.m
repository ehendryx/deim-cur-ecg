function peak_indices = find_peaks(signal_times,signal)
% This function locates local peaks within signals, with parameters
% adjusted heuristically on synthetic data.
%
% Input:    signal_times - a vector containing the timing of each sample in
%                          signal
%           signal - an mxn matrix containing the signals in which to detect
%                    peaks such that each column corresponds to a single signal
%
% Output:   peak_indices - an nx1 cell array containing the corresponding 
%                          indices at which peaks occur in the columns of signal

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the peak threshold to be a multiple of the signal 2-norm scaled by
% the signal length
peak_threshold = 2.5*norm(signal)/sqrt(length(signal));

% Calculate the sampling frequency
dt = diff(signal_times);
mean_dt = mean(dt);
sample_freq = 1/mean_dt;

% Window size in which only one peak can exist
window_size = round(sample_freq/3);

% number of points needed to make up a peak region
min_num_pts = 1;

[row,col]=size(signal);

indices = 1:1:row;

peak_indices = cell(col,1);

% Detect peaks in each column of signal
for counter= 1:col
    
    y = signal(:,counter);
    max_data=[]; % initialize vector of peak indices
    
    % Find those indices greater than the peak threshold
    peak_data_indices = y>peak_threshold;
    
    % Determine where peak regions start and stop
    group_break_indices_up = find(diff(peak_data_indices)>0); % transition into peak region
    group_break_indices_down = find(diff(peak_data_indices)<0); % transition out of peak region
    
    % Require that first peak region starts where the signal is increasing
    if(group_break_indices_down(1) < group_break_indices_up(1))
        group_break_indices_down(1)=[];  % peaks must start on an up mark
    end
    
    % Require that last peak region ends where the signal is decreasing
    if(group_break_indices_up(end) > group_break_indices_down(end))
        group_break_indices_up(end)=[];  % peak must end on a down mark
    end
    
    % Search for peaks in each peak region
    for(group_counter=1:length(group_break_indices_up))
        block_indices = [group_break_indices_up(group_counter):group_break_indices_down(group_counter)];
        data_block_y = y(block_indices);
        data_block_x = indices(block_indices)';
        
        % if the peak region is large enough, identify the peak as the
        % maximimum point in the region
        if(length(data_block_y)>min_num_pts)
            [c,max_index] = max(data_block_y);
            peak_max_location = data_block_x(max_index);
            max_data = [max_data peak_max_location]; % add peak index to list of peaks
        end
    end
    
    peak_indices{counter}=max_data;
    
    
    
    % Eliminate smaller peaks in window regions with multiple peaks
    j = 1;
    while j <= length(peak_indices{counter})
        % define window region to begin at current peak and identify
        % neighboring peaks within this window
        start_time = signal_times(peak_indices{counter}(j));
        end_time_index = min([peak_indices{counter}(j)+window_size-1,length(signal_times)]);
        end_time = signal_times(end_time_index);
        window_peak_indices = find(signal_times(peak_indices{counter})>=start_time & signal_times(peak_indices{counter})<=end_time);
        
        
        if length(window_peak_indices) == 1
            % The jth peak is isolated enough to be kept as is
            j = j+1;
            continue
        else
            % Isolate the highest peak within the window and define that to be the jth peak
            [val,final_peak_index] = max(signal(peak_indices{counter}(window_peak_indices)));
            final_peak_index = j+final_peak_index - 1;
            
            peak_indices{counter}(j) = peak_indices{counter}(final_peak_index);
            peak_indices{counter}(j+1:j+length(window_peak_indices)- 1)=[];
        end
    end
    
    
end
