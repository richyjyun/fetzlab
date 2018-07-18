function [data, end_sec] = nc3chan(start_sec, width_sec, sample_rate, filter_band, output_rate, file_name)
% [data end_sec] = nc3chan(start_sec, width_sec, sample_rate, filter_band, output_rate, file_name)
% Returns a section of analog data from a Neurochip3 data folder.
% Will decimate data if output_rate < channel sample rate (user needs to
% supply appropriate filter_band).
%
% start_sec -- starting time of the data section in seconds
% width_sec -- Width of the data to return.
% sample_rate -- The original sampling rate of data stored in the file.
%   If omitted, a sampling rate of 20000 samples per second is assumed.
%   Use 100 for accelerometer channels.
% filter_band -- [low_cut high_cut] for butterworth filter.  [0 0] for no
%	filter.  May be empty or omitted for no filter.
% output_rate -- Desired sample rate of the resulting data. Maybe omitted,
%   empty, or 0 for no down sampling. If no high_cut is given and
%   output_rate < sample_rate, then high_cut is set to 0.5 * output_rate.
% file_name -- partial or full path name to a neurochip3 .i16 channel file.
%   if omitted, a file dialog box will be displayed to select the .i16 file.
%   A .u16 file is assumed to be and Aux channel.
%
% data -- Returns a column of double precision data. Returns empty on error
%   which includes overreaching the end of file.
% end_sec -- optional, returns the total number of seconds in the file.
%
% Example grabs first 10 seconds of data from a file selected by the user.
% Note that the user must know the sampling rate (20000 in this case).
% >> [data total_seconds_in_file] = nc3chan(0, 10, 20000);
% >> figure; plot(data);
%
% Example grabs first 10 seconds from a file selected by the user, bandpass
% filters [10Hz to 200Hz] and down samples from 20000 to 5000 samples per
% second.
% >> data = nc3chan(0, 10, 20000, [10 200], 5000);

% Default return values

data = [];
end_sec = 0;

% Check parameters

% argin(3) is original sampling rate
inrate = 20000;
if nargin >= 3
    inrate = sample_rate;
end

% argin(4) is filter band.
if (nargin < 4) || (length(filter_band) <= 0)
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

% argin(5) is output sample rate
outrate = inrate;
if nargin >= 5
    outrate = output_rate;
end
if isempty(outrate) || (outrate < 1)
    outrate = inrate;
end

% argin(6) is file name.
fname = '';
if nargin >= 6
    fname = file_name;
end

if (nargin < 2) || (start_sec < 0) || (width_sec <= 0) || (inrate <= 0)
    disp('usage: data = nc3chan(start_sec, width_sec, sample_rate, filter_band, output_rate, file_name)');
    return;
end

% Check if file is a .mat file

if length(fname) < 4
    % Request file from user
    [uifname, uipath] = uigetfile('*.i16', 'Select Neurochip3 .i16 or .u16 Channel file');
    if (uifname(1) == 0) || (uipath(1) == 0)
        return
    end
    fname = fullfile(uipath, uifname);

    % Check for obvious errors.
    if length(fname) < 4
        disp(['Error: ' fname ' is not a Neurochip3 .i16 channel file']);
        return;
    end
end

% Check sample rate and output rate errors.

divide_down = inrate / outrate;
if divide_down ~= floor(divide_down)
    disp('Error: output rate does not divide evenly into sample rate');
    return;
end
chan_samples = width_sec * inrate;
if chan_samples ~= floor(chan_samples)
    disp('Error: time width would give a fractional number of samples at the given sample rate');
    return;
end

data_points = width_sec * outrate;
if (data_points ~= floor(data_points))
    disp('Error: time width would give a fractional number of data points at the given output rate');
    return;
end

% Extract data

adu2units = 0.195;  % From Intan spec sheet.  Intan chip range is 1.25 V per 32768 A/D units with 192 gain.
chtype = 'int16';
if strcmpi(fname(end-3:end), '.u16')
    adu2units = 2.4 / 65536;  % Aux channels return units in Volts
    chtype = 'uint16';
end   
if strcmpi(fname(end-9:end-5), 'Accel')
    adu2units = 1;  % Accel channels are already scaled in milliGrav or Celcius.
end

file_offset = 2 * round(start_sec * inrate);
fid = fopen(fname, 'r');
if (fid ~= -1)
    if nargout > 1
        fseek(fid, 0, 1);
        end_sec = ftell(fid) / (2 * inrate);
    end
    
    if fseek(fid, file_offset, -1) == 0
        % Successfully reached starting sample
        
        [data, count] = fread(fid, [chan_samples, 1], chtype);
        if count ~= chan_samples
            data = [];  % Not enough samples left in file.
        else
            data = adu2units * data;  % convert to uV
            
            % Filter data
            
            b = [1 0 0]; % Default identity filter.
            a = b;

            high = high_cut_Hz;        
            fs2 = 0.5 * inrate;            
%            if (outrate < inrate) && (high == 0)
%                high = 0.5 * outrate;  % Default high cut when none is given but decimation is required.
%            end

            if low_cut_Hz > 0 % Use bandpass filter
                [b, a] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high)] / fs2);             
                data = data - mean(data); % Remove offset to reduce filter edge effects.
                data = filter(b, a, data);
            elseif high > 0 % Use lowpass filter
                [b, a] = butter(2, min(fs2-1, high) / fs2);
                data = filter(b, a, data);
            end

            % Decimate samples if necessary

            if (divide_down > 1)
                data = data(1:divide_down:end);
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