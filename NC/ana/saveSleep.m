function saveSleep(fpath,fname,session_time,sleep_Chn)
% Summary
%   Saves delta oscillations (Delta.mat) and sleep times (Sleep.mat). Saved
%   in fpath.
%
% Inputs
%   fpath - path of the experiment
%   fname - path of the matlab settings file
%   session_time - length of experiment
%   sleep_Chn - channel to check for delta oscillations
%
% RJY 06/22/2018

    % Get delta waves and save
    [Delta,deltafs] = getOscillPower(fname,session_time,sleep_Chn,[1,4]);
    save(fullfile(fpath,'Delta'),'Delta','deltafs');

    % get sleep times and save
    [sleep,sleepfs] = getSleepTimes(fpath,fname,session_time);
    save(fullfile(fpath,'Sleep'),'sleep','sleepfs');
    
end