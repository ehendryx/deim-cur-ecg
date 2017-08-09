function filtered_data = filter_data(signal_times,signal)
% This function applies a basic Butterworth filter

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = butter(1,.005,'high');

filtered_data = filtfilt(b,a,signal);
