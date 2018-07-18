function [sta x n names] = nc3sta_divdown(trig_times, targ_chan_ids, pre_trig_offset_ms, post_trig_offset_ms, bins_per_second, filter_band, file_name)
% [STAs x n names] = nc3sta_divdown(trig_times, targ_chan_ids, pre_trig_offset_ms, post_trig_offset_ms, bins_per_second, filter_band, file_name)
% Calculates a spike triggered average between the trigger times and the
% target channels.
%
% trig_times -- event times given in seconds from the beginning of the file.
% targ_chan_id -- list of channel ID numbers.  Use [0:15] for first 16 channels.
% pre_trig_offset_ms -- analysis pre-trigger time offset in milliseconds.
%   This is usualy  negative so that the average encompasses time zero.
% post_trig_offset_ms -- analysis post-trigger time offset in milliseconds.
%   This is usualy positive so that the average encompasses time zero.
% bins_per_second -- Number of bins per second in the resulting average.
% filter_band -- [low_cut high_cut] for butterworth filter.  [0 0] for no
%	filter.  May be omitted for no filter.
% file_name -- partial or full path name to the neurochip3 .mat settings
%   file in the data folder. if omitted, a file dialog box will be
%   displayed to select the .mat file.
%
% sta -- Returns the spike triggered average y values in microvolts.
% x -- Returns the x values (in seconds) for each bin.
% n -- Actual number of trig_times used for each average.
% names -- Returns cell array of channel names.
%
% Example:
% >> [STAs x n] = nc3sta_divdown(trigs, [0:15], -50, 100, 10000, [500 2500]);
% >> figure;
% >> plot(x, STAs);  % plot all spike triggered averages stacked.
% >> figure;
% >> plot(x * 1000, STAs(:,1));  % plot the first average in milliseconds.
% >> title(['STA of ' names{1} ' (n = ' num2str(n(1)) ')']);
% >> xlabel('ms');
% >> ylabel('uV');

% Default return value

sta = [];
x = [];
n = [];
names = {};

% Check parameters

if nargin < 5
    outrate = 5000;
else
    outrate = bins_per_second;
end
if (outrate < 100)
    outrate = 100;
    disp('Warning: Bins_per_second set to 100');
end
fs2 = outrate / 2;

if nargin < 6
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

if nargin < 7
    fname = '';
else
    fname = file_name;
end

nchans = length(targ_chan_ids);
ntrigs = length(trig_times);
if (ntrigs <= 0) || (nchans <= 0) || (pre_trig_offset_ms > post_trig_offset_ms) || (high_cut_Hz > fs2)
    disp('usage: [STAs x n names] = nc3sta(trig_times, targ_chan_ids, pre_trig_offset_ms, post_trig_offset_ms, bins_per_second, filter_band, file_name)');
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
    [b a] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high_cut_Hz)] / fs2);
elseif high_cut_Hz > low_cut_Hz % Use lowpass filter
    [b a] = butter(2, min(fs2-1, high_cut_Hz) / fs2);
end

% Check sample rates, output rate, and filter band for errors.

rates = zeros(nchans, 1);
divide_down = ones(nchans, 1);
preoffset_samples = zeros(nchans,1);
sweep_samples = zeros(nchans,1);
for ichan = 1:nchans
    id = targ_chan_ids(ichan);
    if (id >= start_index) && ((id + 1 - start_index) <= length(p.channel_rate)) && (p.channel_rate(id + 1 - start_index) > 0)
        names{ichan} = p.channel_names{id + 1 - start_index};
        rates(ichan) = p.channel_rate(id + 1 - start_index);
        divide_down(ichan) = rates(ichan) / outrate;
        if divide_down(ichan) ~= floor(divide_down(ichan))
            disp(['Error: bins_per_second does not divide evenly into sample rate for channel ' num2str(id)]);
            return;
        end
        preoffset_samples(ichan) = (pre_trig_offset_ms * rates(ichan) / 1000);
        sweep_samples(ichan) = ((post_trig_offset_ms - pre_trig_offset_ms) * rates(ichan) / 1000);
        if (preoffset_samples(ichan) ~= floor(preoffset_samples(ichan))) || (sweep_samples(ichan) ~= floor(sweep_samples(ichan)))
            disp(['Error: Pre-trigger and post-trigger offsets cannot be divided evenly by sample rate for channel ' num2str(id)]);
            return;
        end
    else
        disp(['Error: channel id ' num2str(id) ' does not exist in ' fname]);
        return;
    end
end

preoffset_bins = pre_trig_offset_ms * outrate / 1000;
sweep_bins = (post_trig_offset_ms - pre_trig_offset_ms) * outrate / 1000;
if (preoffset_bins ~= floor(preoffset_bins)) || (sweep_bins ~= floor(sweep_bins))
    disp('Error: Pre-trigger and post-trigger offsets cannot be divided evenly by bins_per_second');
    return;
end

% Compile aveages

sta = zeros(sweep_bins, nchans);
n = zeros(nchans,1);
adu2uv = 1250000 / (32768 * 196);  % Intan chip range is 1.25 V per 32768 A/D units with 196 gain.

for ichan = 1:nchans
    if rates(ichan) > 0
        ave = zeros(sweep_bins, 1);
        file_offsets = 2 * round(trig_times * rates(ichan) + preoffset_samples(ichan));
        nsweeps = 0;
        id = targ_chan_ids(ichan);
        cname = [fpath fprefix '_Chan' num2str(id, '%.02d') '.i16'];
        fid = fopen(cname, 'r');
        start_trig = find(file_offsets >= 0, 1);
        if (fid ~= -1) && ~isempty(start_trig)
            for itrig = start_trig:ntrigs
                if fseek(fid, file_offsets(itrig), -1) ~= 0
                    break; % Seek failed, we must be outside of file limits.
                end

                [dat, count] = fread(fid, [sweep_samples(ichan), 1], 'int16');
                if count ~= sweep_samples(ichan)
                    break; % Could not read the entire sweep.
                end

                if (divide_down(ichan) > 1)
                    % Decimate samples down to bins.
                    dat = mean(reshape(dat, divide_down(ichan), sweep_bins), 1)';
                end

                if (low_cut_Hz > 0)
                    % Remove offset to reduce filter edge effects.
                    dat = dat - mean(dat);
                end

                if (high_cut_Hz > 0)
                    dat = filter(b, a, dat);
                end

                ave = ave + dat;
                nsweeps = nsweeps + 1;
            end

            sta(:, ichan) = ave * adu2uv / nsweeps;  % convert to uV
            fclose(fid);
        end
        n(ichan) = nsweeps;
    end
end

% X axis in seconds

x = ((0:sweep_bins-1) + preoffset_bins) / outrate;

% example figure

if 1
    figure;
    for ichan = 1:nchans
        subplot(nchans, 1, ichan)
        plot(x * 1000, sta(:,ichan))
        ylabel('uV');
        xlabel(['STA ' names{ichan} ' (n = ' num2str(n(ichan)) ')']);
        if ichan == 1
            title(fname);
        end
    end
end