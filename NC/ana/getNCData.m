function [fpath,fname,Channels,fs,session_time] = getNCata(path,day)
% Summary
%   Returns relevant meta data extracted from the NC data
%
% Inputs
%   path - path of all NC data
%   day - specific experiment
%
% Outputs
%   fpath - path of the experiment
%   fname - name of the matlab file
%   Channels - recorded channels
%   fs - sampling rate
%   session_time - length of recording
%
% RJY 06/22/2018

    fpath = fullfile(path,day);
    fname = fullfile(path,day,[day,'.mat']);

    % Error checking
    if(~exist(path))
        error('The data path does not exist');
    elseif(~exist(fpath))
        error('The experiment does not exist');s
    elseif(~exist(fname))
        error('The Matlab setting file does not exist')
    end

    Events = nc3events(fname);
    Channels = find(Events.channel_rate(1:32) > 0);
    fs = Events.channel_rate(Channels);
    
    % if channels are recorded at different sampling rates, send to
    % keyboard for user to resolve
    if(length(unique(fs)) > 1)
        warning('Channels have different sampling rates')
        keyboard;
    end
    
    fs = unique(fs);

    session_time = Events.session_time;

end