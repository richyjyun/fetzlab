function [sleep,fs] = getSleepTimes(fpath,fname,session_time)
% Summary
%   Determines when animal is sleeping using delta power and accelerometer
% data
%
% Inputs
%   fpath           path of the NC data
%   fname           NC settings file name
%   session_time    length of NC session in seconds
% 
% Ouputs
%   sleep           time determined to be sleep (in samples)
%   fs              sampling rate
%
% RJY 06/22/2018

    % static downsample frequency
    fs = 100;

    %% Load Delta power
    load(fullfile(fpath,'Delta.mat'));

    % get movinng average and normalize
    d = movmean(Delta,100*1000); 
    d = d-min(d);

    % set arbitrary threshold
    thresh = std(d);

    % find times above threshold
    highd = d >= thresh; 

    %% Load z-axis accelerometer data
    Accel = nc3data(38, 0, session_time-mod(session_time,fs), fs, [], fname);
    Accel = interp(Accel,10);
    
    % get moving average
    a = movmean(Accel,100*1000);
    thresh = 300; % arbitrary, always going to be the same with accelerometer

    % find times below threshold
    lowa = a <= thresh;

    %% Define sleep time
    def = min([length(highd),length(lowa)]);
    sleep = highd(1:def) & lowa(1:def);
    
    fs = 1000;

end



