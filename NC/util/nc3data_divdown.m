function [data, names] = nc3data_divdown(chan_ids, start_sec, width_sec, output_rate, filter_band, file_name)
% [data names] = nc3data(chan_ids, start_sec, width_sec, output_rate, filter_band, file_name)
% Returns a section of analog data from a Neurochip3 data folder.
%
% chan_ids -- list of channel ID numbers. Use [0:15] for pre v320 files,
% or [1..32] for v320+ files
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
% >> [data, names] = nc3data_divdown([0:15], 0, 10, 5000, [500 2500]);
% >> plot(data(:,1));  % plot the first channel

% Default return value

data = [];
names = {};

% Check parameters

if nargin < 4
    outrate = 5000;
else
    outrate = output_rate;
end
if (outrate < 100)
    outrate = 100;
    disp('Warning: output rate set to 100 samples per second');
end
fs2 = outrate / 2;

if nargin < 5
    low_cut_Hz = 0;
    high_cut_Hz = 0;
else
    if length(filter_band) < 2
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
if (nargin < 3) || (nchans <= 0) || (start_sec < 0) || (width_sec <= 0) || (high_cut_Hz > fs2)
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

% Calculate filter coefficients

b = [1 0 0]; % Default identity filter.
a = b;

if low_cut_Hz > 0 % Use bandpass filter
    [b, a] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high_cut_Hz)] / fs2);
elseif high_cut_Hz > low_cut_Hz % Use lowpass filter
    [b, a] = butter(2, min(fs2-1, high_cut_Hz) / fs2);
end

% Check sample rates, output rate, and filter band for errors.

rates = zeros(nchans, 1);
divide_down = ones(nchans, 1);
chan_samples = zeros(nchans,1);
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

adu2uv = 0.195;  % From Intan spec sheet.  Intan chip range is 1.25 V per 32768 A/D units with ~192 gain.

for ichan = 1:nchans
    if (rates(ichan) > 0)
        file_offset = 2 * round(start_sec * rates(ichan));
        id = chan_ids(ichan);
        cname = [fpath fprefix '_Chan' num2str(id,'%.02d') '.i16'];
        fid = fopen(cname, 'r');
        if (fid ~= -1)
            if fseek(fid, file_offset, -1) ~= 0
                break; % Seek failed, we must be outside of file limits.
            end

            [dat, count] = fread(fid, [chan_samples(ichan), 1], 'int16');
            if count ~= chan_samples(ichan)
                break; % Could not read the entire sweep.
            end

            if (divide_down(ichan) > 1)
                % Decimate samples down to bins.
                dat = mean(reshape(dat, divide_down(ichan), data_points), 1)';
            end

            if (low_cut_Hz > 0)
                % Remove offset to reduce filter edge effects.
                dat = dat - mean(dat);
            end

            if (high_cut_Hz > 0)
                dat = filter(b, a, dat);
            end

            if isempty(data)
                data = zeros(data_points, nchans);
            end
            data(:, ichan) = dat * adu2uv;  % convert to uV
            fclose(fid);
        end
    end
end

% example figure

if 0
    figure;
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