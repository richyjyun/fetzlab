function [data, end_sec] = nc3chandata(start_sec, width_sec, sample_rate, bad_sections, margin_samples, data_range, file_name)
% [data end_sec] = nc3chandata(start_sec, width_sec, sample_rate, bad_sections, margin_samples, data_range, file_name)
% Returns a section of analog data from a Neurochip3 data folder.
% Interpolates over bad sections of data.  Bad sections can be explicitly
% given or defined implicitly by defining a range of good data samples.
% Bad ranges are extended backwards and forwards by margin_samples before
% interpolation.
%
% start_sec -- starting time of the data section in seconds
% width_sec -- Width of the data to return.
% sample_rate -- The original sampling rate of data stored in the file.
%   If omitted, a sampling rate of 20000 samples per second is assumed.
% bad_sections -- a list of missing or bad secttions usually from the
%   nc3events() functions. This is a 2xN matrix of time intervals where
%   events.missing_sections(1,n) is the first unusable time (in seconds).
%   events.missing_sections(2,n) is the next usuable time (in seconds). 
%   if empty or omitted, then there are no missing sections of data.
% margin_samples -- A margin for bad data removal given in data samples.
%   0 if omitted.  All bad sections are extended by margin_samples on
%   both ends.
% data_range -- Acceptable data range in uV, data outside this range is
%   flagged as bad and will be removed by linear interpolation between
%   a margin before and after the the bad data.  all data ok if omitted.
% file_name -- partial or full path name to a neurochip3 .i16 channel file.
%   if omitted, a file dialog box will be displayed to select the .i16 file.
%
% data -- Returns a column of double precision data. Returns empty on error
%   which includes overreaching the end of file.
% end_sec -- optional, returns the total number of seconds in the file.
%
% Example grabs first 10 seconds of data from a file selected by the user.
% All missing sections and values outside of +/-100 uV are interpolated
% over. Note that the user must know the sampling rate (20000 in this case)
%
% >> events = nc3events(experiment_filename);
% >> [data, total_seconds_in_file] = nc3chandata(0, 10, 20000, events.missing_sections, 1, [-100 100], channel_filename);
% >> figure; plot(data);

% Default return values

data = [];
end_sec = 0;

% Check parameters

% argin(3) is original sampling rate
inrate = 20000;
if nargin >= 3
    inrate = sample_rate;
end

% argin(4) is bad time sections
if (nargin < 4)
    bad_sections = [];
end

% argin(5) is margin samples
if nargin < 5
    margin_samples = 0;
end
margin_samples = floor(max(0, margin_samples));

% argin(6) is the acceptable data range
if nargin < 6
    data_range = [];
end

% argin(7) is file name.
fname = '';
if nargin >= 7
    fname = file_name;
end

if (nargin < 2) || (start_sec < 0) || (width_sec <= 0) || (inrate <= 0)
    disp('usage: nc3chandata(start_sec, width_sec, sample_rate, bad_sections, margin_samples, data_range, file_name)');
    return;
end

% Check if file is a .mat file

if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.i16') == 0)
    % Request file from user
    [uifname, uipath] = uigetfile('*.i16', 'Select Neurochip3 .i16 Channel file');
    if (uifname(1) == 0) || (uipath(1) == 0)
        return
    end
    fname = fullfile(uipath, uifname);

    % Check for obvious errors.
    if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.i16') == 0)
        disp(['Error: ' fname ' is not a Neurochip3 .i16 channel file']);
        return;
    end
end

% Check sample rate errors.

chan_samples = width_sec * inrate;
if chan_samples ~= floor(chan_samples)
    disp('Error: time width would give a fractional number of samples at the given sample rate');
    return;
end

% Extract data

adu2uv = 0.195;  % From Intan spec sheet.  Intan chip range is 1.25 V per 32768 A/D units with ~192 gain.

first_sample = round(start_sec * inrate);
file_offset = 2 * first_sample;
fid = fopen(fname, 'r');
if (fid ~= -1)
    if nargout > 1
        fseek(fid, 0, 1);
        end_sec = ftell(fid) / (2 * inrate);
    end
    
    if fseek(fid, file_offset, -1) == 0
        % Successfully reached starting sample
        
        [data, count] = fread(fid, [chan_samples, 1], 'int16');
        if count ~= chan_samples
            data = [];  % Not enough samples left in file.
        else
            data = adu2uv * data;  % convert to uV
            
            % Create flags for marking bad intervals
            
            bad_flags = zeros(chan_samples, 1, 'int8');
            
            % Mark out samples in the bad_sections list
            
            [r,c] = size(bad_sections);
            if (r == 2) && (c > 0)
                bad_samples = floor(bad_sections * inrate) - first_sample + 1;
                bad_index = find((bad_samples(2,:) >= 1) & (bad_samples(1,:) <= chan_samples));
                n = length(bad_index);
                for i = 1:n
                    s1 = max(1, bad_samples(1,bad_index(i)) - margin_samples);
                    s2 = min(chan_samples, bad_samples(2,bad_index(i)) + margin_samples);
                    bad_flags(s1:s2) = 1;
                end
            end
            
            % Mark out samples outside the range of good data

            if length(data_range) >= 2
                bad_index = find((data < data_range(1)) | (data > data_range(2)));
                n = length(bad_index);
                for i = 1:n
                    s1 = max(1, bad_index(i) - margin_samples);
                    s2 = min(chan_samples, bad_index(i) + margin_samples);
                    bad_flags(s1:s2) = 1;
                end
            end
            
            % Linear interpolation over bad intervals
            
            bad_start = find((bad_flags(1:end-1) == 0) & (bad_flags(2:end) == 1));
            bad_stop = find((bad_flags(1:end-1) == 1) & (bad_flags(2:end) == 0));
            
            if bad_flags(1) == 1 % Special case a leading interval
                if isempty(bad_stop)
                    data(1:end) = mean(data); % entire data section is bad
                else
                    data(1:bad_stop(1)) = data(bad_stop(1));
                    bad_stop(1) = [];  % Remove first stop.
                end
            end
            
            if bad_flags(end) == 1 && ~isempty(bad_start) % Special case a trailing interval
                data(bad_start(end):end) = data(bad_start(end)); 
                bad_start(end) = []; % Remove last start
            end
            
            n = min(length(bad_start), length(bad_stop));
            for i = 1:n
                s1 = bad_start(i);  % Bad interval start sample index.
                s2 = bad_stop(i);   % stop sample index
                width = s2-s1;      % Interval width.
                if (width > 1)
                    val1 = data(s1); % Starting value
                    val2 = data(s2); % Stopping value
                    step_size = (val2 - val1) / width;
                    data(s1:s2) = (step_size * (0:width)) + val1;  % Linear interpolation between val1 .. val2
                end
            end
        end
    end
    fclose(fid);
end

% Example figure

if 0
    if isempty(data)
        disp('Data returned empty -- could not read requested section of data');
    else
        figure;
        plot(data);
        ylabel('uV');
        xlabel(fname);
    end
end