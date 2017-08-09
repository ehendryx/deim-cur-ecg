function data_matrix = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude)
% This function generates a matrix containing individual beats detected 
% from a synthetic ECG time series. The time series is generated given a set of
% input parameters for the desired feature characteristics within the
% series.
%
% Input:    N - number of R-R intervals to generate
%           HR - 1x2 array containing the desired signal heart rate and the
%                  maximum amount of increase allowed
%           Rpeak - shift parameter from the center of each beat time interval
%           *amp - 1x2 array containing the corresponding desired feature
%                  amplitude parameter and the maximum amount of variability allowed
%           *std - 1x2 array containing the corresponding desired feature
%                  amplitude parameter and the maximum amount of variability allowed
%           *shift - 1x2 array containing the corresponding desired feature
%                  shift parameter relative to the R peak and the maximum
%                  amount of variability allowed; note that T_shift and P_shift are scaled by 
%                  the beat time in the code below; also note that the resulting
%                  shifts in defining feature location (except for RRshift) must be nonnegative
%           baseline - 0 indicates no baseline wandering, 1 adds sinusoidal baseline wandering
%           baseline_mag - a parameter to determine the scale of baseline wandering
%           series_magnitude - a parameter for scaling the magnitude of the generated time series ,file_description     
%   See end of code for some example parameter values.
%   Note:Parameter order for *amp, *std, and *shift: [value, max % change allowed]
%
% Output:   data_matrix - an Nx125 matrix with each row containing an individual 
%                         beats from the generated synthetic ECG time series

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

N = N+1; % accounts for fact that individual beats will be generated starting with P-wave and ending with T-wave.

% Sampling frequency
fs = 240; 

% Add variability to heart rate
HR = HR(1)*(1+rand(N,1)*HR(2));

% Determine the beat centers from the HR of each beat 
t0 = zeros(1,N+1);
t0_index = zeros(1,N+1);

t0(1) = 0;
t0_index(1) = 1;
beat_time = 60/HR(1); % Convert to seconds per beat
t = linspace(t0(1),t0(1) + beat_time,round(beat_time*fs)); % initialize full ECG series time samples
for i = 2:N
    beat_time = 60/HR(i); % Convert to seconds per beat
    t0(i) = t(end); % Begin new beat where old beat ended
    t0_index(i) = length(t);
    new_t = linspace(t0(i),t0(i) + beat_time,round(beat_time*fs));
    t = [t(1:end-1), new_t]; % Add new beat time samples to full ECG series 
end
t0(N+1) = t(end);
t0_index(N+1) = length(t);


% Generate baseline shift
if baseline == 0
    baseline = zeros(1,length(t));
else
    baseline = baseline_mag*(.03*sin(2*pi/(4*beat_time)*t) + .05*sin(2*pi/(10*beat_time)*t) + .03*sin(2*pi/beat_time*t));
end


% Generate individual beats
beat = cell(N,1);

j = 1;
while j <= N 
    time_indices = t0_index(j):t0_index(j+1); % global time indices
    local_time = t(time_indices); % times for individual beat
    [flag, beat{j,1}] = gaussian_generator_v2(local_time,t,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd);
    
    if ~isempty(beat{j,1})
        j = j+1;
    end
    
    % Make sure beat features were generated in the appropriate order
    if flag == 1
        return
    end
end

% Sum individual beats to make ECG series
series = beat{1,1};

for k = 2:N
    series = series + beat{k,1};
end

% Add baseline and scale by series_magnitude
series = series + baseline;

series = series_magnitude*series;


% Save beats from R peak to R peak

peak_indices = find_peaks(t(:),abs(series(:)));

% Make sure the desired number of beats were detected
if numel(peak_indices{1}) ~= N
    display('Incorrect number of peaks detected. Run again, and consider using different parameters.')
    return
end


% Construct Beat Matrix

d = 125; % number of samples stored from each beat through interpolation
data_matrix = NaN(length(peak_indices{1})-1,d);

for i = 1:length(peak_indices{1})-1
    interp_times  = linspace(t(peak_indices{1}(i)),t(peak_indices{1}(i+1)),d);
    interp_data = interp1(t(peak_indices{1}(i):peak_indices{1}(i+1)),series(peak_indices{1}(i):peak_indices{1}(i+1)),interp_times,'linear');
    data_matrix(i,:) = interp_data;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example input values:
% HR = [120, HR_percent_change];
%
% Rpeak = 0; % This entry is a temporal shift from center of beat.
% Ramp = [1, amp_percent_change];
% Rstd = [1e-2, std_percent_change];
% 
% Pshift = [.25, shift_percent_change]; % This entry is a fraction of the total beat the total beat time to shift from the R peak. (Example: %Pshift(1)*beat_time*(1+rand*percent_change)
% Pamp = [.1, amp_percent_change];
% Pstd = [0.02, std_percent_change];
% 
% Tshift = [.3, shift_percent_change]; %.25; %This entry is a fraction of the total beat the total beat time to shift from the R peak. (Example: %Tshift(1)*beat_time*(1+rand*percent_change)
% Tamp = [.15, amp_percent_change];
% Tstd = [.03, std_percent_change];
% 
% Qshift = [Rstd(1), shift_percent_change]; %Rstd(1)*mvnrnd(1,.01);
% Qamp = [.2, amp_percent_change];
% Qstd = [1e-2, std_percent_change];
% 
% Sshift = [Rstd(1), shift_percent_change]; %Rstd(1)*mvnrnd(1,.01);
% Samp = [.25, amp_percent_change];
% Sstd = [1e-2, std_percent_change];
% 
% RRshift = [0, shift_percent_change];
% RRamp = [0, amp_percent_change]; % no R' peak
% RRstd = [.8*Rstd(1), std_percent_change];
% 
% series_magnitude = 500/400;

