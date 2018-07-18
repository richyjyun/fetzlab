function nc3filter(low_cut_Hz, high_cut_Hz, sample_rate, down_sample, file_name)
% nc3filter(low_cut_Hz, high_cut_Hz, sample_rate, down_sample, file_name)
% Filters and down_samples a Neurochip3 data file.
% Parameters down_sample and file_name are optional.
%
% [low_cut_Hz, high_cut_Hz] -- the filter bandwidth.
% sample_rate -- The channel's sample rate which can be found in the .mat
%   settings file or the .txt info file (20000, 10000, or 5000)
% down_sample -- Output rate = sample_rate / down_sample.  Defaults to 1.
% file_name -- partial or full path name to a neurchip3 analog channel file
%   This file must have a .i16 file extension.
%
% This function creates an output file '<file_pathname>_<low>_<high>_<rate>.i16
% which stores the best estimate of the filter used in the Neurochip3
% window discriminator for the given bandwidth and sample rate.
% The output file is truncated to a whole number of seconds.
% If no input file is given, a dialog box will be displayed to ask for one.

if nargin < 4
    divdown = 1;
else
    divdown = down_sample;
end

if nargin < 5
    fname = '';
else
    fname = file_name;
end

fs2 = sample_rate / 2;

% Check parameters.
if (nargin < 3) || (divdown < 1) || (low_cut_Hz < 0) || (high_cut_Hz < low_cut_Hz) || (high_cut_Hz > fs2) || (fs2 < 2500) || (fs2 > 10000)
    disp('Usage: nc3filter(low_cut_Hz, high_cut_Hz, sample_rate, down_sample, file_name)');
    disp('Parameters down_sample and file_name are optional.');
    return;
end

% Check if file is a .i16 file

if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.i16') == 0)
    % Request file from user
    [uifname, uipath] = uigetfile('*.i16', 'Select Neurchip3 .i16 data file to Convert to filter');
    if (uifname(1) == 0) || (uipath(1) == 0)
        return
    end
    fname = fullfile(uipath, uifname);

    % Check for obvious errors.
    if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.i16') == 0)
        disp(['Error: ' fname ' is not a Neurochip3 .i16 data file']);
        return;
    end
end

% Check if we can open the data file.

fid = fopen(fname, 'r');
if fid < 0
    disp(['Error: Could not open ' fname]);
    return;
end

% Check if we can open the output file.

outrate = sample_rate / divdown;
outname = [fname(1:end-4) '_f_' num2str(low_cut_Hz) '_' num2str(high_cut_Hz) '_' num2str(outrate) '.i16'];
outfid = fopen(outname, 'w');
if outfid < 0
    disp(['Error: Could not write to ' outname]);
    fclose(fid);
    return;
end

% Calculate filter coefficients

b = [1 0 0]; % Default identity filter.
a = b;

if low_cut_Hz > 0 % Use bandpass filter
    [b, a] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high_cut_Hz)] / fs2);
elseif high_cut_Hz > low_cut_Hz % Use lowpass filter
    [b, a] = butter(2, min(fs2-1, high_cut_Hz) / fs2);
end

% Quantize filter coefficients to 24 bits to match Neurochip3 event filters.
% Comment this out if floating point precision should be maintained.

%a = round(a * (2^24))./(2^24);
%b = round(b * (2^24))./(2^24);

% Initialize filter state

z = zeros(max(length(a),length(b)) - 1, 1);

% Scan through file 1 second at a time

seconds = 0;
while 1
    % Read input samples in whole seconds
    [data, count] = fread(fid, sample_rate, 'int16');
    if count ~= sample_rate
        break;
    end
    
    % Filter data
    if low_cut_Hz < high_cut_Hz
        [data, z] = filter(b, a, data, z);
    end
    
    % Down sample data
    if (divdown > 1)
        data = mean(reshape(data, divdown, outrate), 1);
    end
    
    % Write output samples
    count = fwrite(outfid, data, 'int16');
    if (count ~= outrate)
        disp(['Error: output file is incomplete: ' outname]);
        break;
    end
    seconds = seconds + 1;
end

% Close files and clean up.

disp([num2str(seconds) ' seconds of data written to: ' outname]);
fclose(fid);
fclose(outfid);
