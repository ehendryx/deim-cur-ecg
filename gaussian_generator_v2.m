function [flag, beat] = gaussian_generator_v2(t,total_time,Rpeak,Ramp,Rstd,Pshift,Pamp,Pstd,Qshift,Qamp,Qstd,Sshift,Samp,Sstd,Tshift,Tamp,Tstd,RRamp,RRshift,RRstd)
% This function generates individual synthetic ECG beats to be added
% together in forming a series of beats. Beats are constructed with
% features defined relative to the R peak.
%
% Input:    t - a vector containing the desired time samples over which the beat should occur
%           total_time - time samples for the full series, allowing Gaussians
%                     to decay naturally rather than abruptly truncating at
%                     the beat's edges
%           Rpeak - shift parameter from the center of the time interval
%                   given in t
%           *amp - 1x2 array containing the corresponding desired feature
%                  amplitude parameter and the maximum amount of variability allowed
%           *std - 1x2 array containing the corresponding desired feature
%                  amplitude parameter and the maximum amount of variability allowed
%           *shift - 1x2 array containing the corresponding desired feature
%                  shift parameter relative to the R peak and the maximum
%                  amount of variability allowed; note that T_shift and P_shift are scaled by 
%                  the beat time in the code below; also note that the resulting
%                  shifts in defining feature location (except for RRshift) must be nonnegative
% Parameter order for *amp, *std, and *shift: [value, max % change allowed]

% Output:   flag - a value of 1 indicates that input parameters resulted in
%                   a negative feature shift for input into the beat-defining
%                   Gaussians, which can result in features placed out of
%                   order
%           beat - the resulting individual beat placed in the context of
%                  total_time

% This code is under a 3-Clause BSD License.
% Copyright 2017, E. Hendryx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag = 0;

beat_time = t(end)-t(1);


% Define parameters for the beat-defining Gaussians based on the input
% values
Rpeak = Rpeak+t(round(length(t)/2));
Ramp = Ramp(1)*(1+(rand*2*Ramp(2) - Ramp(2))); % The window of allowed variability is centered around 0.
Rstd = Rstd(1)*(1+(rand*2*Rstd(2) - Rstd(2)));

Pshift = (Pshift(1)*beat_time)*(1+(rand*2*Pshift(2) - Pshift(2))/2);
Pamp = Pamp(1)*(1+(rand*2*Pamp(2) - Pamp(2)));
Pstd = Pstd(1)*(1+(rand*2*Pstd(2) - Pstd(2)));

Tshift = (Tshift(1)*beat_time)*(1+(rand*2*Tshift(2) - Tshift(2))/2);
Tamp = Tamp(1)*(1+(rand*2*Tamp(2) - Tamp(2)));
Tstd = Tstd(1)*(1+(rand*2*Tstd(2) - Tstd(2)));

Qshift = Qshift(1)*(1+(rand*2*Qshift(2) - Qshift(2))/2);
Qamp = Qamp(1)*(1+(rand*2*Qamp(2) - Qamp(2)));
Qstd = Qstd(1)*(1+(rand*2*Qstd(2) - Qstd(2)));

Sshift = Sshift(1)*(1+(rand*2*Sshift(2) - Sshift(2))/2);
Samp = Samp(1)*(1+(rand*2*Samp(2) - Samp(2)));
Sstd = Sstd(1)*(1+(rand*2*Sstd(2) - Sstd(2)));

RRshift = RRshift(1)*(1+(rand*2*RRshift(2) - RRshift(2))/2);
RRamp = RRamp(1)*(1+(rand*2*RRamp(2) - RRamp(2)));
RRstd = RRstd(1)*(1+(rand*2*RRstd(2) - RRstd(2)));


% Make sure features will appear on the correct side of the R peak given
% the signs used in the beat-defining Gaussians
if Pshift < 0 || Qshift < 0 || Sshift < 0 || Tshift < 0
    display('All shifts, except for RRshift, must be nonnegative')
    display(['Pshift = ' num2str(Pshift)])
    display(['Qshift = ' num2str(Qshift)])
    display(['Sshift = ' num2str(Sshift)])
    display(['Tshift = ' num2str(Tshift)])
    flag = 1;
    return
end

% Define individual features based on the above parameters
R = Ramp*exp(-((total_time - Rpeak)/Rstd).^2);

Q = -Qamp*exp(-((total_time - (min(Rpeak,Rpeak+RRshift)-Qshift))/Qstd).^2);

S = -Samp*exp(-((total_time - (max(Rpeak,Rpeak+RRshift)+Sshift))/Sstd).^2);

P = Pamp*exp(-((total_time - (min(Rpeak,Rpeak+RRshift)-Pshift))/Pstd).^2);

T = Tamp*exp(-((total_time - (max(Rpeak,Rpeak+RRshift)+Tshift))/Tstd).^2);

RR = RRamp*exp(-((total_time - (Rpeak+RRshift))/RRstd).^2);

% sum the features to form a full beat
beat = sparse(P+Q+R+S+T+RR);
