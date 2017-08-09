function data_matrix = trial_matrix_generator(HR_percent_change,shift_percent_change,amp_percent_change,std_percent_change,noise_percentage)
% This function generates a data matrix containing individual 6000 beats from 12
% different synthetic ECG series (500 beats from each series). Each row of
% the data matrix contains a beat and is Z-normalized.

% Input:    HR_percent_change - amount of variation allowed in the
%                               synthetic heart rate given as a fraction
%           shift_percent_change - amount of variation allowed in the
%                               synthetic feature positions given as a
%                               fraction
%           amp_percent_change - amount of variation allowed in the
%                               synthetic feature amplitudes given as a
%                               fraction
%           std_percent_change - amount of variation allowed in the
%                               synthetic heart rate given as a fraction
%           noise_percentage - amount of noise added to the data matrix 
%                               given as a fraction 
%
% Output:   data_matrix - a 6000x125 matrix containing 500 from each of 12
%                       synthetically generated ECG time series; each row
%                       contains an individual Z-normalized beat

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of beats to generate in each class type
N = 500;

% Do not include baseline wandering in signal generation
baseline = 0;
baseline_mag = 0;

% Define control heart rate and allowed variability
HR = [120, HR_percent_change];

% Initialize data matrix
data_matrix = [];

%% Generate 12 different classes of heart beats %%
% Parameter order: [value, percent_change]

% See gaussian_generator_v2.m to understand how these parameters are used
% in defining beats.

% Class 1:
% Standard features present (no R' peak)
class = 1;

Rpeak = 0;
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.3, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.2, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.25, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 2:
% No S-wave, and an R' peak is present
class = class + 1;

Rpeak = 0; 
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.3, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change]; 
Qamp = [.25, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change]; 
Samp = [0, amp_percent_change]; % no S wave
Sstd = [1e-2, std_percent_change];

RRshift = [2.25*Rstd(1), shift_percent_change];
RRamp = [.5, amp_percent_change];
RRstd = [.6*Rstd(1), std_percent_change];

series_magnitude = 600/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 3:
% No P-wave
class = class + 1;

Rpeak = 0;
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [0, amp_percent_change]; % no P-wave
Pstd = [0.02, std_percent_change];

Tshift = [.3, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.2, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.25, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 700/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 4:
% Standard features are present and R' is also present.
class = class + 1;

Rpeak = 0;
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.3, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.25, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.2, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [2.25*Rstd(1), shift_percent_change];
RRamp = [.7, amp_percent_change];
RRstd = [.6*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 5:
% Standard features present (no R' peak), but T-wave is premature looking
class = class + 1;

Rpeak = 0;
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.15, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.25, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.4, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 650/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 6:
% Premature T wave with no R peak and more pronounced S wave (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [0, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.15, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.2, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.5, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 7:
% Slightly premature, inverted T wave with no R peak and wider, inverted S wave, making it look kind of like there's a bi-modal T wave (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [0, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.04, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.25, shift_percent_change];
Tamp = [-.05, amp_percent_change];
Tstd = [.05, std_percent_change];

Qshift = [1.1*Rstd(1), shift_percent_change];
Qamp = [.35, amp_percent_change];
Qstd = [0.01, std_percent_change];

Sshift = [1.3*Rstd(1), shift_percent_change];
Samp = [-.02, amp_percent_change];
Sstd = [0.07, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 600/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 8:
% Small T wave with wider S wave (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [.4, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.04, amp_percent_change];
Pstd = [0.03, std_percent_change];

Tshift = [.25, shift_percent_change];
Tamp = [.005, amp_percent_change];
Tstd = [.05, std_percent_change];

Qshift = [1.1*Rstd(1), shift_percent_change];
Qamp = [.05, amp_percent_change];
Qstd = [0.01, std_percent_change];

Sshift = [3*Rstd(1), shift_percent_change];
Samp = [.02, amp_percent_change];
Sstd = [0.08, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 700/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 9:
% Inverted T wave with wider S wave and no Q wave (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [.4, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.02, amp_percent_change];
Pstd = [0.03, std_percent_change];

Tshift = [.28, shift_percent_change];
Tamp = [-.06, amp_percent_change];
Tstd = [.04, std_percent_change];

Qshift = [1.1*Rstd(1), shift_percent_change];
Qamp = [0, amp_percent_change];
Qstd = [0.01, std_percent_change];

Sshift = [4*Rstd(1), shift_percent_change];
Samp = [.02, amp_percent_change];
Sstd = [0.06, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 700/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 10:
% Small R peak with prominent S wave and no Q wave (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [.2, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.25, shift_percent_change];
Pamp = [.02, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.3, shift_percent_change];
Tamp = [.03, amp_percent_change];
Tstd = [.035, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [0, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.3, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 11:
% % Small R peak with prominent S wave and no Q wave, as well as overlapping T and P waves (no R' peak)
class = class + 1;

Rpeak = 0;
Ramp = [.2, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.45, shift_percent_change];
Pamp = [.02, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.4, shift_percent_change];
Tamp = [.03, amp_percent_change];
Tstd = [.035, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [0, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.3, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [0, shift_percent_change];
RRamp = [0, amp_percent_change]; % no R' peak
RRstd = [.8*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)

%% Class 12:
% Standard features are present and R' is also present; T and P waves are overlapping
class = class + 1;

Rpeak = 0;
Ramp = [1, amp_percent_change];
Rstd = [1e-2, std_percent_change];

Pshift = [.4, shift_percent_change]; 
Pamp = [.1, amp_percent_change];
Pstd = [0.02, std_percent_change];

Tshift = [.4, shift_percent_change];
Tamp = [.15, amp_percent_change];
Tstd = [.03, std_percent_change];

Qshift = [Rstd(1), shift_percent_change];
Qamp = [.25, amp_percent_change];
Qstd = [1e-2, std_percent_change];

Sshift = [Rstd(1), shift_percent_change];
Samp = [.2, amp_percent_change];
Sstd = [1e-2, std_percent_change];

RRshift = [2.25*Rstd(1), shift_percent_change];
RRamp = [.7, amp_percent_change];
RRstd = [.6*Rstd(1), std_percent_change];

series_magnitude = 500/400;

new_data = ECG_series_generator_v2(N,HR,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd,baseline,baseline_mag,series_magnitude); 

data_matrix = [data_matrix; new_data];

figure
plot(max_pos)
hold on
plot(max_neg,'k')
title(['Class ' num2str(class)],'fontsize',14)


%%% Add Noise %%%
if noise_percentage ~= 0
    noise_matrix = randn(size(data_matrix));
    noise_matrix = noise_matrix/norm(noise_matrix)*norm(data_matrix);
    
    data_matrix = data_matrix + noise_percentage*noise_matrix;
end


%%% Normalize %%%
data_matrix = data_matrix_beat_normalization(data_matrix);





