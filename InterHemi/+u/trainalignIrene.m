function error_code = trainalignIrene(fileprefix, trig_index, train_mV_threshold)
% Requires the following files
% <fileprefix>.i16 -- The saved data stream from the training computer.
% <fileprefix>.bin -- The data stream from daqdiscrim/daqbimanual programs
% <fileprefix>.cfg -- The config file from daqdiscrim/daqbimanual
% trig_index -- 1 or 2.  Use first or second trigger channel for alignment.
% train_mv_threshold -- Optional trigger threshold in mV. 3000 = default.
%
% Saves <fileprefix.f32> as an aligned version of <fileprefix>.i16.
% The .i16 data is upsampled and aligned to match the data in the .bin file.
% Alignment is done by matching a pulse train on the seventh channel
% in the .i16 file to one of the trigger channels in the .bin file.
% Upsampling is done by linear interpolation between adjacent samples
% for the first 6 channels, and earliest neighbor for the last 3 channels.
% The output file always contains 9 channels written as 'single'
% precision as a [9 x N] matrix.  N should match the number of samples in the
% .bin file, and the first and last sample may be repeated as padding to
% cover any time shift between the train and guger data files.
%
% error_code:
% 0 = Output file created (although there may be some displayed warnings)
% 1 = Could not find the .i16 training file
% 2 = Could not find the .bin data file.
% 3 = Could not find the .cfg config file.
% 4 = Could not open the output file for writing.
% 5 = trig_index was not 1 or 2, or .cfg file did not have the trigger.
% 7 = Fewer than 10 triggers were found in the .bin file.
% 8 = Fewer than 10 triggers were found in the .i16 file.
% 9 = Triggers could not be matched reliably.

fnname = 'trainalignIrene';

train_threshold = 3000;
if nargin >= 3
    train_threshold = train_mV_threshold;
end

if (trig_index < 1) || (trig_index > 2)
    disp(['Error: ' fnname '. trig_index must be 1 or 2']);
    error_code = 5;
    return
end

infile = [fileprefix '.i16'];
datfile = [fileprefix '.bin'];
cfgfile = [fileprefix '.cfg'];
outfile = [fileprefix '.f32'];

[infile, error_code] = MakeTempI16(fileprefix, infile); % added by ARB 12/17/14

infid = fopen(infile, 'r');
if (infid < 0)
    disp(['Error: ' fnname ' could not find training data file' infile]);
    error_code = 1;
    return
end
datfid = fopen(datfile, 'r');
if (datfid < 0)
    fclose(infid);
    disp(['Error: ' fnname ' could not find daq data file ' datfile]);
    error_code = 2;
    return
end
cfgfid = fopen(cfgfile, 'r');
if (cfgfid < 0)
    fclose(infid);
    fclose(datfid);
    disp(['Error: ' fnname ' could not find daq config file ' cfgfile]);
    error_code = 3;
    return
end
fclose(cfgfid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read config file

load(cfgfile,'-mat');
nchans = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
if isfield(UI, 'wd_channel')
    nchans = nchans + 1; % Window discriminator indicates daqdiscrim file.
else
    nchans = nchans + 2; % Assume daqbinmanual was used.
end
chnm = cell(nchans,1); % name of channels
chgu = zeros(nchans,1); % gUSBamp that channel was recorded on
ii = 0;
for iga = 1:length(UI.ga_trigger)
    ind = find(UI.ch_enabled(iga,:));
    chnm(ii+1:ii+length(ind)) = UI.ch_name(iga,ind)';
    chgu(ii+1:ii+length(ind)) = iga*ones(length(ind),1);
    ii = ii + length(ind);
    if UI.ga_trigger(iga), chnm{ii+1} = ['trig ' num2str(iga)]; chgu(ii+1) = iga; ii = ii + 1; end
end
if isfield(UI, 'wd_channel')
    chnm{end} = 'Discrim';  % One extra channel in daqdiscrim files.
else
    chnm{end-1} = 'Behave1'; % Two extra chanels in daqbimanual files.
    chnm{end} = 'Behave2';
end
itrigch = strmatch('trig',chnm);
fs = UI.samprate;

if length(itrigch) < trig_index
    fclose(infid);
    fclose(datfid);
    disp(['Error: ' fnname ' could not find specified trigger channel in ' cfgfile]);
    error_code = 5;
    return    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read triggers from .bin file

trig_offset = (itrigch(trig_index)-1) * 4;
fseek(datfid, trig_offset, -1);
skip_bytes = (nchans - 1) * 4;
[samples total_samples] = fread(datfid, '*single', skip_bytes);
fclose(datfid);
trig1 = find(samples > 1);       % Find triggers
intervals1 = diff([0; trig1]);   % Calcualte intervals
trig1 = (trig1(intervals1 > 10) - 1) * 1000 / fs; % Debounce trigger and remove doublets. Time base of zero.
intervals1 = diff([0; trig1]); % Intervals in milliseconds.
trig1_mean = mean(intervals1(2:end));
disp(['Mean data trigger interval = ' num2str(trig1_mean) ' ms (N=' num2str(length(trig1)) ')']);
if (length(trig1) < 10)
    disp(['Error: ' fnname ' found less than 10 triggers in ' datfile]);
    fclose(infid);
    error_code = 7;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read triggers from .i16 file

trig_offset = (7 - 1) * 2;  % offset to 7th channel
fseek(infid, trig_offset, -1);
skip_bytes = (9 - 1) * 2;   % 9 signals, 8 of which need to be skipped.
[samples train_samples] = fread(infid, 'int16=>single', skip_bytes);
index = find(samples > train_threshold); % look for 3 consecutive samples above our threshold
trig2 = index((index(1:end-2)+1 == index(2:end-1)) & (index(2:end-1)+1 == index(3:end)));
intervals2 = diff([0; trig2]);
trig2 = trig2(intervals2 > 2)-1;% Debounce trigger and remove doublets. Time base of zero.
intervals2 = diff([0; trig2]); % intervals in milli seconds.
trig2_mean = mean(intervals2(2:end));
disp(['Mean train trigger interval = ' num2str(trig2_mean) ' ms (N=' num2str(length(trig2)) ')']);
if (length(trig2) < 10)
    disp(['Error: ' fnname ' found less than 10 triggers in ' infile]);
    fclose(infid);
    error_code = 8;
    return;
end

% Look for a good interval match at the beginning of the file

n1 = length(trig1);
n2 = length(trig2);

match1 = 0;
match2 = 0;
match_score = inf;
matchlimit = min([10 (n1-2) (n2-2)]);
for i1=1:matchlimit
    for i2 = 1:matchlimit
       % score = i1 + i2 + sum(abs(intervals1(i1:i1+2) - intervals2(i2:i2+2)));
        score = abs(trig1(i1) - trig2(i2));
        if score < match_score
            match_score = score;
            match1 = i1;
            match2 = i2;
        end
    end
end

offset = trig1(match1) - trig2(match2);
disp(['Guger trigger ' num2str(match1) ' at ' num2str(floor(trig1(match1))) ...
    ' ms matches Train trigger ' num2str(match2) ' at ' num2str((trig2(match2))) ...
    ' ms.  Offset = ' num2str((offset)) ' ms.']);

% Look for a good time match in the middle of the file.

n = min([n2 - 7, n1 - 7, floor(n1/2)]);
endmatch_score = inf;
endmatch1 = 0;
endmatch2 = 0;
for i1 = n:n+5
    [score, index] = min(abs(trig2 + offset - trig1(i1)));
    score = score + sum(abs(intervals1(i1:i1+2) - intervals2(index:index+2)));
    if score < endmatch_score
        endmatch_score = score;
        endmatch1 = i1;
        endmatch2 = index;
    end
end

endoffset = trig1(endmatch1) - trig2(endmatch2);
disp(['Guger trigger ' num2str(endmatch1) ' at ' num2str(floor(trig1(endmatch1))) ...
    ' ms matches Train trigger ' num2str(endmatch2) ' at ' num2str((trig2(endmatch2))) ...
    ' ms.  Offset = ' num2str((endoffset)) ' ms.']);

clock_skew = double(offset - endoffset) / double(trig1(endmatch1) - trig1(match1));  % Revised sample rate of train file.
disp(['Adjusting for clock skew of ' num2str(clock_skew * 60 * 60 * 1000) ' ms per hour of recording.']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the output file

outfid = fopen(outfile, 'w');
if (outfid < 0)
    fclose(infid);
    disp(['Error: ' fnname ' could not save output file ' outfile]);
    error_code = 4;
    return
end

% Create output file 100000 samples at a time.

for isamp = 0:100000:total_samples-1
    % Find fractional sample times in train file adjusting for clock skew.
    xtimes = (((isamp:min(total_samples-1, isamp+99999)) * 1000 / fs) - trig1(match1)) * (1 + clock_skew) + trig2(match2);
    xindex1 = floor(xtimes);   % First nearest sample index.
    xfrac2 = xtimes - xindex1; % Fractional overlap with the second sample.
    xfrac1 = 1 - xfrac2;       % Fractional overlap with first sample.
    % Samples occuring before beginning or after end of file will
    % be adjusted to use the first or last sample in the file.
    if (xindex1(1) < 1)
        xindex1(xindex1 < 0) = 0;
    end
    if (xindex1(end) > train_samples-1)
        xindex1(xindex1 > train_samples-1) = train_samples-1;
    end
    fseek(infid, xindex1(1) * 18, -1); % offset to first sample used
    n = xindex1(end) - xindex1(1) + 2; % number of samples referenced.
    samples = fread(infid, [9 n], 'int16');
    xindex1 = xindex1 - xindex1(1) + 1; % Convert from time to block sample index
    xindex2 = xindex1 + 1;
    n = size(samples,2);
    if xindex2(end) > n
        xindex2(xindex2 > n) = n; % adjust for end of file.
    end
    outbuf = zeros(9, length(xindex1));
    for i=1:9
        if i<=6
            outbuf(i,:) = (xfrac1 .* samples(i,xindex1)) + (xfrac2 .* samples(i,xindex2)); % interpolate analog signals
        else
            outbuf(i,:) = samples(i, xindex1); % Earliest neighbor for digital signals
        end
    end
    fwrite(outfid, outbuf, 'single');
end

fclose(infid);
fclose(outfid);
error_code = 0;

%%%%%%%%%%%%%%%%%%%%%%
% Verify output file by reading in the aligned trigger channel in the .f32 file.

disp(['Verifying triggers in ' outfile]);
infid = fopen(outfile, 'r');
trig_offset = (7 - 1) * 4;  % offset to 7th channel
fseek(infid, trig_offset, -1);
skip_bytes = (9 - 1) * 4;   % 9 signals, 8 of which need to be skipped.
[samples train_samples] = fread(infid, 'single=>single', skip_bytes);
index = find(samples > train_threshold); % look for 3 consecutive samples above our threshold
trig2 = index((index(1:end-2)+1 == index(2:end-1)) & (index(2:end-1)+1 == index(3:end)));
intervals2 = diff([0; trig2]);
trig2 = (trig2(intervals2 > 2)-1)* 1000 / fs; % Debounce trigger and remove doublets. Time base of zero.
fclose(infid);

% Compare aligned triggers with guger triggers

n = 0;
sum_offsets = 0;
n_offsets = 0;
for i1 = 1:n1
    [score index] = min(abs(trig2 - trig1(i1)));
    if score > 2
        disp(['Guger trigger at ' num2str(trig1(i1)) ' ms failed to match a Train trigger within 2ms.  Nearest is offset by ' num2str(trig2(index) - trig1(i1)) 'ms.']);
        n = n + 1;
    else
        n_offsets = n_offsets + 1;
        sum_offsets = sum_offsets + score;
    end
end

% Display final stats.

disp(['Total time in .bin = ' num2str(total_samples * 1000 / fs) ' ms.']);
disp(['Samples in .bin = ' num2str(total_samples)]);
disp(['Samples in .f32 = ' num2str(train_samples)]);
disp(['Average matching trigger jitter = ' num2str(sum_offsets / n_offsets) ' ms.']);
if n == 0
    disp('Success! All triggers match within 2ms.');
else
    disp(['Warning: '  num2str(n) ' Guger triggers had no match within 2ms.']);
end
if (n > 10) || (total_samples ~= train_samples)
    disp('!!! Verify found issues.  Please check results !!!');
    error_code = 9;
end

delete tmp.i16 % in this format so that a real i16 isnt accidentally deleted ARB 12/17/14

function [tmpinfile, error_code] = MakeTempI16(fileprefix, infile)
% This function makes a new, corrected i16 file that will be deleted after trainalign.m runs
%
% andrew 17 dec 2014

fs = 1000; % this is always the same for i16, as far as I know

tmpinfile = 'tmp.i16'; % temporary i16 file

fid = fopen(infile, 'r');
if (fid < 0)
    disp(['Error: ' fnname ' could not find training data file' infile]);
    error_code = 1;
    return
end
fdata = fread(fid, [9, inf], 'int16=>single');
fclose(fid);

% added by arb 12/15/14
if str2double(fileprefix(end-10:end-3))<20141208 % these have a relatively systematic offset between binary and analogue data
    disp('Data file is from before Dec 8th 2014, an assumed shift within training file will be updated. Know that training program timings such as box presentation times will be imperfectly aligned to neural/behavioral data with a jitter +/- ~5ms')
    outbuf = fdata(9,:); % digital data that needs to be shifted forward by 100ms  
    shiftsamples = round(fs/10); % 100ms worth of shift
    outbuf = [repmat(outbuf(1), 1, shiftsamples), outbuf(1:end-shiftsamples)]; % pad with first value 
    fdata(9,:) = outbuf;
elseif str2double(fileprefix(end-10:end-3))>=20141208
    disp('Data file is from after Dec 8th 2014. Shift within training file will be precisely corrected.')
    
    targetonscreen = bitand(uint16(fdata(9,:)), 2); % this is correctly timed
    boxid = floor(fdata(9,:) / 1024); % this is used for timing reference
    
    ev_boxid = findOnsetsAndOffsets(boxid>0); % calculate periods for reference timing
    ev_tos = findOnsetsAndOffsets(targetonscreen>0);
    
    [~, ~, ds_on] = FindClosestBefore(ev_tos(:,1), ev_boxid(:,1)); % find offsets in time
    [~, ~, ds_off] = FindClosestBefore(ev_tos(:,2), ev_boxid(:,2));
    
    ds = [ds_on, ds_off];
    ds = ds(2:end-1,:); % clip the first and last pairs
    ev_tos = ev_tos(2:end-1,:);

    offsets = interp1(reshape(ev_tos', numel(ev_tos), 1), reshape(ds', numel(ds),1), 1:length(boxid)); % interpolate offsets in time

    inds = (1:length(boxid)) - round(offsets); % these represent the indices in box id that make up the corrected indices from first ev_ts to last ev_ts 
    inds(isnan(inds)) = [];
    
    bitstoshift = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 16]; % shift all but the target on screen bits

    bindat = dec2bin(fdata(9,:),16);
    bindat(ev_tos(1):ev_tos(end), bitstoshift) = bindat(inds, bitstoshift);
    bindat(1:ev_tos(1)-1, bitstoshift) = repmat(bindat(inds(1), bitstoshift), ev_tos(1)-1, 1); % buffer front
    bindat(ev_tos(end)+1:end, bitstoshift) = repmat(bindat(inds(end), bitstoshift), size(bindat,1) - ev_tos(end), 1); % buffer end
    fdata(9,:) = bin2dec(bindat); % update binary data    
end

fid = fopen(tmpinfile, 'w'); % write file
fwrite(fid, single(fdata), 'int16');
fclose(fid);

fid = fopen(tmpinfile, 'r'); % checks
fdata = fread(fid, [9, inf], 'int16=>single');
fclose(fid);
fid = fopen(infile, 'r'); % checks
fdata2 = fread(fid, [9, inf], 'int16=>single');
fclose(fid);

error_code = 0;

function [Y, inds, ds] = FindClosestBefore(x, y)
% [Y inds ds] = FindClosest(x, y)
% 
% This function finds the value in vector Y that is closest and less than
% to every value in X. Values that are equal are discarded.
%
% Returns
%   Y: a vector where Y(i) is taken from Y and minimizes ds = X(i)-Y(i) for ds>0
%   inds: indices in Y that match above
%   ds: a vector == x-Y
%
% arb 11/10/14

y = sort(y);

inds = ones(size(x));

for i = 1:length(x)
   
    tmp = find(x(i)-y>0, 1, 'last');
    
    if ~isempty(tmp), inds(i) = tmp; end
    
end

ds = x-y(inds);

Y = y(inds); % the values in sorted Y

Y(ds<1) = NaN;
inds(ds<1) = NaN;
ds(ds<1) = NaN;

function [OnOffs] = findOnsetsAndOffsets(boolVec)
% function [OnOffs] = findOnsetsAndOffsets(boolVec)
%
% Returns list of aligned start and stops of chunks of 1's in a vector of
% 1's and 0's.
%
% INPUTS
%  boolVec - vector of 1's and 0's 
%
% OUTPUTS
%  startEnds - Nx2 list of indices of the first and last 1's for the N
%              contiguous blocks of 1's.
%

% created 6/16/11 eln
% update 7/5/11 arb line 1
% update 7/6/11 eln line 26

boolVec = boolVec(:)';

starts = find(diff(boolVec)==1);
ends = find(diff(boolVec)==-1);

% if starts out going correct speed, add 1 to starts
if boolVec(1)
  starts = [0 starts];
end

% if finishes going correct speed, add final value to ends
if boolVec(end)
  ends = [ends length(boolVec)];
end

OnOffs = [starts(:)+1 ends(:)];
