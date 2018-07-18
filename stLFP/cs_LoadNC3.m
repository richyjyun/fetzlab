clear; close all;

fpath = 'Y:\~NeuroWest\Spanky\Neurochip\S_20180605_01\';
chn = 27; fs = 20000;

% [data, names, session_time] = nc3data(chn, 0, 100, 20000, [300,3000], path);

index = strfind(fpath, '\');
fprefix = fpath(index(end-1)+1:end-1);
fname = [fpath,fprefix,'_Chan' num2str(chn,'%.02d') '.i16'];

% initialize sorting parameters
init_time = 3600; % time to use 
[Data, end_sec] = nc3chan(0, 3600, fs, [], fs, fname);

[p,spk,noise] = getSortingParams(Data,fs);

% sort and find spike times
step = 3600; % look at an hour at a time
times = step:step:end_sec;
times = [times,end_sec];

for t = times
    [Data, end_sec] = nc3chan(t-step, t, fs, [], fs, fname);

    
    
end
