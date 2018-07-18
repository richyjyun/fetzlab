function nc3combine(low_cut_Hz, high_cut_Hz, output_rate, file_name)
% nc3combine(low_cut_Hz, high_cut_Hz, output_rate, file_name)
% Filters and possibly decimates Neurochip3 data to an interleaved
% output file.  file_name of the .mat file is optional, but if omitted,
% a file dialog box will be displayed to select the .mat file.
%
% [low_cut_Hz, high_cut_Hz] -- the filter bandwidth.
% output_rate -- The output sample rate for the interleaved output file.
% file_name -- partial or full path name to the neurchip3 .mat settings file in the data folder.
%
% This function creates an output file '<file_pathname>_<low>_<high>_<num_chans>_<output_rate>.i16
% which stores a filtered and interleaved version of the all the analog channel data.
% If low_cut_Hz is 0, a high pass filter is used (otherwise a band pass filter is used)
% If high_cut_Hz is 0, no filter is used unless the output_rate is lower than a
% channel's sample rate, then a high_cut is set at half the output_rate to reduce
% aliasing in that channel's decimated data.


if nargin < 4
    fname = '';
else
    fname = file_name;
end

fs2 = output_rate / 2;

% Check parameters.
if (nargin < 3) || (low_cut_Hz < 0) || (high_cut_Hz < low_cut_Hz) || (high_cut_Hz > fs2) || (output_rate > 20000)
    disp('Usage: nc3filter(low_cut_Hz, high_cut_Hz, output_rate, file_name)');
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

% Open upto 16 input file.
% Make sure output_rate divides eveningly into all sample rates.

infids = zeros(32,1);
rates = zeros(32,1);
num_files = 0;
for ifile = 1:32
    cname = [fpath fprefix '_Chan' num2str(ifile-1,'%.02d') '.i16'];
    fid = fopen(cname, 'r');
    if fid ~= -1
        divdown = rates(ifile) / output_rate;
        if divdown ~= floor(divdown)
            disp(['Warning: divide down of ' num2str(divdown) ' is not supported. Skipping ' p.channel_names{ifile}]);
            fclose(fid);
            continue
        end
        num_files = num_files + 1;
        infids(num_files) = fid;
        rates(num_files) = p.channel_rate(ifile);
        decimation{num_files} = divdown;
        
        b{num_files} = [1 0 0]; % Default identity filter.
        a{num_files} = b{num_files};
        z{num_files} = zeros(max(length(a),length(b)) - 1, 1); % Initialize filter state

        high = high_cut_Hz;
        fs2 = 0.5 * rates(num_files);
        if (output_rate < rates(num_files)) && (high == 0)
            high = 0.5 * output_rate;  % Default high cut when none is given but decimation is required.
        end

        if low_cut_Hz > 0 % Use bandpass filter
            [b{num_files}, a{num_files}] = butter(1, [min(fs2-1, low_cut_Hz) min(fs2-1, high)] / fs2);
        elseif high > 0 % Use lowpass filter
            [b{num_files}, a{num_files}] = butter(2, min(fs2-1, high) / fs2);
        end

    end
end

if num_files == 0
    disp(['Error: no channels to combine for ' fname]);
    return;
end

% Check if we can open the output file.

outname = [fname(1:end-4) '_' num2str(low_cut_Hz) '_' num2str(high_cut_Hz) '_' num2str(num_files) '_' num2str(output_rate) '.i16'];
outfid = fopen(outname, 'w');
if outfid < 0
    disp(['Error: Could not write to ' outname]);
    for ifile = 1:num_files
        fclose(infids(ifile));
    end
    return;
end

% Allocate output buffer

outbuf = zeros(num_files, output_rate);
outbuf_size = output_rate * num_files;

% Scan through file 1 second at a time

seconds = 0;
while 1
    % Read input samples in whole seconds
    for ifile = 1:num_files
        [data, count] = fread(infids(ifile), rates(ifile), 'int16');
        if count ~= rates(ifile)
            count = 0;
            disp(num2str(count));
            break;
        end
        
        % Filter data 
        if low_cut_Hz < high_cut_Hz
            [data, z{ifile}] = filter(b{ifile}, a{ifile}, data, z{ifile});
        end

        % Decimate data if needed.
        if (decimation{ifile} > 1)
            data = data(1:decimation{ifile}:end);
        end

        % fill in each row of the output buffer.
        outbuf(ifile, :) = data;
    end
    if count == 0
        break;
    end

       
    % Write interleaved output samples.
    count = fwrite(outfid, outbuf, 'int16');
    if (count ~= outbuf_size)
        disp(['Error: output file is incomplete. Out of disk space?: ' outname]);
        break;
    end
    seconds = seconds + 1;
end

% Close files and clean up.

fclose(outfid);
for ifile = 1:num_files
    fclose(infids(ifile));
end

disp([num2str(seconds) ' seconds of data written to: ' outname]);
