function [Oscill,fs] = getOscillPower(fname,session_time,chn,bp)
% Summary
% Returns oscillation power of defined frequency band in the specific channel
% 
% Inputs
% fname - NC settings file name
% session_time - length of NC session in seconds
% chn - NC channel to load
% bp - frequency band (e.g. [15,25])
% 
% Outputs
% Oscill - oscillations
% fs - sampling rate
%
% RJY 06/22/2018

    % static downsample frequency
    fs = 1000;

    %% Load data 
    Oscill = [];

    data = nc3data(chn, 0, floor(session_time), fs, bp, fname);
    Oscill = abs(hilbert(data));

end

