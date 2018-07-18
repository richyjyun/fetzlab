function [data, names, session_time] = nc3data(chan_ids, start_sec, width_sec, output_rate, filter_band, file_name)
% [data names session_time] = nc3data(chan_ids, start_sec, width_sec, output_rate, filter_band, file_name)
% Returns a section of analog data from a Neurochip3 data folder.
% Will decimate data if output_rate < channel sample rate (user needs to
% supply appropriate filter_band)
% Will pad data requested past the end of file with zeros.
% session_time is optional and returns the total number of seconds in the file.
% Biophysical channels are returned in microVolts, Aux channels in Volts,
% accelerometer data (at 100 samples per second) is in units of 0.001 gravity.
%
% chan_ids -- list of channel ID numbers.  Use [0..31] for pre v320 files,
%   or [1..32] for v320+ files. Aux channels are [32..34] for pre v32 files,
%   or [33..35] for v320+ files. Accelerometer channels follow Aux channels
%   with AccelX, AccelY, AccelZ, AccelT, AccelM.  
% start_sec -- starting time of the data section in seconds
% width_sec -- Width of the data to return.
% output_rate -- Desired sample rate of the resulting data.
% filter_band -- [low_cut high_cut] for butterworth filter.  [0 0] for no
%	filter.  May be omitted for no filter.
% file_name -- partial or full path name to the neurochip3 .mat settings
%   file in the data folder. if omitted, a file dialog box will be
%   displayed to select the .mat file.
%
% data -- Returns the data in columns.  Returns empty on error.
% names -- Returns cell array of channel names.
%
% Example grabs first 10 seconds from a file selected by the user
% >> [data names] = nc3data([0:15], 0, 10, 5000, [500 2000]);
% >> plot(data(:,1));  % plot the first channel

% Default return value

data = [];
names = {};
if nargout > 2
    session_time = 0;
end

% Check parameters

if nargin < 4
    outrate = 20000;
else
    outrate = output_rate;
end
if (outrate < 100)
    outrate = 100;
    disp('Warning: output rate set to 100 samples per second');
end

if (nargin < 5) || (length(filter_band) <= 0)
    low_cut_Hz = 0;
    high_cut_Hz = 0;
else
    if length(filter_band) <= 1
        low_cut_Hz = 0;
        high_cut_Hz = filter_band(1);
    else
        low_cut_Hz = filter_band(1);
        high_cut_Hz = filter_band(2);
    end
end

if nargin < 6
    fname = '';
else
    fname = file_name;
end

nchans = length(chan_ids);
if (nargin < 3) || (nchans <= 0) || (start_sec < 0) || (width_sec <= 0) || (high_cut_Hz > 0.5 * output_rate)
    disp('usage: [data names] = nc3data(chan_ids, start_sec, width_sec, output_rate, filter_band, file_name)');
    return;
end

% Check if file is a .mat file

if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.mat') == 0)
    % Request file from user
    [uifname, uipath] = uigetfile('*.mat', 'Select Neurchip3 Settings file');
    if (uifname(1) == 0) || (uipath(1) == 0)
        return
    end
    fname = fullfile(uipath, uifname);

    % Check for obvious errors.
    if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.mat') == 0)
        disp(['Error: ' fname ' is not a Neurochip3 settings file']);
        return;
    end
end

% Load neurochip3 parameters.

try
    index = strfind(fname, '\');
    if isempty(index)
        fpath = '';
        fprefix = fname(1:end-4);
    else
        fpath = fname(1:index(end));
        fprefix = fname(index(end)+1:end-4);
    end
    loadstruct = load(fname);
    p = loadstruct.p;
catch
    disp(['Error: ' fname ' is not a Neurochip3 settings file']);
    return;
end

start_index = 0;
if p.version >= 320
    start_index = 1; % Version 320 files: first channel index is 1.
end

% Check sample rates, output rate, and filter band for errors.

rates = zeros(nchans, 1);
divide_down = ones(nchans, 1);
chan_samples = zeros(nchans,1);
names = cell(nchans, 1);

for ichan = 1:nchans
    id = chan_ids(ichan);
    if (id >= start_index) && ((id + 1 - start_index) <= length(p.channel_rate)) && (p.channel_rate(id + 1 - start_index) > 0)
        names{ichan} = p.channel_names{id + 1 - start_index};
        rates(ichan) = p.channel_rate(id + 1 - start_index);
        divide_down(ichan) = rates(ichan) / outrate;
        if divide_down(ichan) ~= floor(divide_down(ichan))
            disp(['Error: output rate does not divide evenly into sample rate for channel ' num2str(id)]);
            return;
        end
        chan_samples(ichan) = width_sec * rates(ichan);
        if chan_samples(ichan) ~= floor(chan_samples(ichan))
            disp(['Error: data time width cannot be divided evenly by sample rate for channel ' num2str(id)]);
            return;
        end
    else
        disp(['Error: channel id ' num2str(id) ' does not exist in ' fname]);
        return;
    end
end

data_points = width_sec * outrate;
if (data_points ~= floor(data_points))
    disp('Error: Pre-trigger and post-trigger offsets cannot be divided evenly by bins_per_second');
    return;
end

% extract data

data = zeros(data_points, nchans);  % Initialize output array.
file_time = 0; % Number of seconds in the first file.

for ichan = 1:nchans
    if (rates(ichan) > 0)
        file_offset = 2 * round(start_sec * rates(ichan));
        id = chan_ids(ichan);
        adu2units = 1; % Default no units conversion.
        idcheck = id - start_index;
        cname = '';
        ctype = 'int16';
        if (idcheck <= 31)
            adu2units = 0.195;  % Adu to uV. From Intan spec sheet. Intan chip range is 1.25 V per 32768 A/D units with ~192 gain.
            cname = [fpath fprefix '_Chan' num2str(id,'%.02d') '.i16'];
        elseif (idcheck <= 34)
            adu2units = 2.4 / 65536;  % Adu to V. Aux channels are 2.4 volts at max 16-bit unsigned value
            cname = [fpath fprefix '_Chan' num2str(id,'%.02d') '.u16'];
            ctype = 'uint16';
        elseif (idcheck == 35)
            cname = [fpath fprefix '_AccelX.i16'];
        elseif (idcheck == 36)
            cname = [fpath fprefix '_AccelY.i16'];
        elseif (idcheck == 37)
            cname = [fpath fprefix '_AccelZ.i16'];
        elseif (idcheck == 38)
            cname = [fpath fprefix '_AccelT.i16'];
        elseif (idcheck == 39)
            cname = [fpath fprefix '_AccelM.i16'];
        end
        fid = fopen(cname, 'r');
        if (fid ~= -1)
            if file_time == 0
                fseek(fid, 0, 1);
                file_time = ftell(fid) / (2 * rates(ichan));
            end

            if fseek(fid, file_offset, -1) ~= 0
                data = [];
                break; % Seek failed, we must be outside of file limits.
            end

            [dat, count] = fread(fid, [chan_samples(ichan), 1], ctype);
            if (count < chan_samples(ichan))      
                dat(count+1:chan_samples(ichan)) = 0; % pad with zeros to correct length
            end
            
            % Filter channel data.

            b = [1 0 0]; % Default identity filter.
            a = b;

            high = high_cut_Hz;        
            fs2 = 0.5 * rates(ichan);            
%            if (outrate < rates(ichan)) && (high == 0)
%                high = 0.5 * outrate;  % Default high cut when none is given but decimation is required.
%            end
            
            if low_cut_Hz > 0 % Use bandpass filter
                [b, a] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high)] / fs2);
            elseif high > 0 % Use lowpass filter
                [b, a] = butter(2, min(fs2-1, high) / fs2);
            end

            if (low_cut_Hz > 0)
                % Remove offset to reduce filter edge effects.
                dat = dat - mean(dat);
            end

            if (high > 0)
                dat = filter(b, a, dat);
            end

            % Decimate samples if necessary
            
            if (divide_down(ichan) > 1)
                dat = dat(1:divide_down(ichan):end);
            end

            data(1:length(dat), ichan) = dat * adu2units;  % convert to proper units
            fclose(fid);
        end
    end
end

if nargout > 2
    session_time = file_time;
end

% example figure

if 0
    figure; %#ok<*UNRCH>
    for ichan = 1:nchans
        subplot(nchans, 1, ichan)
        plot(data(:,ichan))
        ylabel('uV');
        xlabel(names{ichan});
        if ichan == 1
            title(fname);
        end
    end
end