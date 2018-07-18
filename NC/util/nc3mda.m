function nc3mda(filename)
% Function nc3mda(filename)
%   filename -- name of the settings file inside the neurochip 3 data folder.
%       If the file name is omitted, a file dialog box is opened.
% Convert neurchip3 data files to MDA format.  MDA files are stored in
% folder c:\data\mda\

datapath = 'c:\data\mda\';  % path to store converted mda files

if nargin < 1
    fname = '';
else
    fname = filename;
end

% Check if file is a .mat file

if (length(fname) < 4) || (strcmpi(fname(end-3:end), '.mat') == 0)
    % Request file from user
    [uifname, uipath] = uigetfile('*.mat', 'Select Neurchip3 Settings file to Convert to MDA');
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

% Setup base MDA object.

experimentName = [fprefix '_mda'];
if isdir(datapath) == 0
    mkdir(datapath);
end
mdadir = [datapath 'Data for ' experimentName '\'];
mkdir(mdadir);

Version = 2;
Objects = []; % Create Root experiment object
Objects(1).PropertyNames{1} = 'Name';
Objects(1).PropertyValues{1} = experimentName;
Objects(1).PropertyNames{2} = 'Animal';
Objects(1).PropertyValues{2} = experimentName(1:strfind(experimentName, '_')-1);
Objects(1).PropertyNames{3} = 'Date';
Objects(1).PropertyValues{3} = date;
Objects(1).PropertyNames{4} = 'DataFilenameExt';
Objects(1).PropertyValues{4} = '';
Objects(1).PropertyNames{5} = 'fileName';
Objects(1).PropertyValues{5} = experimentName;
Objects(1).PropertyNames{6} = 'fileNameTemplate';
Objects(1).PropertyValues{6} = '<Name>';
Objects(1).PropertyNames{7} = 'DataPath';
Objects(1).PropertyValues{7} = ['~\Data for ' experimentName];
Objects(1).PropertyNames{8} = 'DataPathTemplate';
Objects(1).PropertyValues{8} = '~\Data for <FullName>';
Objects(1).PropertyNames{9} = 'FullName';
Objects(1).PropertyValues{9} = experimentName;
Objects(1).PropertyNames{10} = 'FullNameTemplate';
Objects(1).PropertyValues{10} = '<Name>';

Objects(1).Parent = 1;
Objects(1).Children = [];
Objects(1).Class = 'experiment';
Objects(1).SubClass = '';
Objects(1).Changed = 1;
Objects(1).Version = 2;
Objects(1).ID = 2;
id = 1;  % Current object ID number.

% Convert analog channels.

experimentTime = 0;
maxTime = 0;
zeroHeader = zeros(1, 1024);

for ifile = 0:40
    if ifile > length(p.channel_rate)
        break;
    end
    if ifile == 0
        % Treat binary channel as an analog at 5000 samples per second.
        cname = [fpath fprefix '_Digi00.u16'];
        rate = 5000;
        name = 'Digital';
    elseif ifile == 36
        cname = [fpath fprefix '_AccelX.i16'];
        rate = 100;
        name = 'AccelX';
    elseif ifile == 37
        cname = [fpath fprefix '_AccelY.i16'];
        rate = 100;
        name = 'AccelY';
    elseif ifile == 38
        cname = [fpath fprefix '_AccelZ.i16'];
        rate = 100;
        name = 'AccelZ';
    elseif ifile == 39
        cname = [fpath fprefix '_AccelM.i16'];
        rate = 100;
        name = 'AccelM';
    elseif ifile == 40
        cname = [fpath fprefix '_AccelT.i16'];
        rate = 100;
        name = 'AccelT';
    else
        % Do next analog channel
        if p.version < 320
            cname = [fpath fprefix '_Chan' num2str(ifile-1,'%.02d') '.i16']; % for nc310
        elseif ifile <= 32
            cname = [fpath fprefix '_Chan' num2str(ifile,'%.02d') '.i16']; % for nc320 normal Intan channels
        else
            cname = [fpath fprefix '_Chan' num2str(ifile,'%.02d') '.u16']; % nc320 AUX channels
        end
        rate = p.channel_rate(ifile);
        name = p.channel_names{ifile};
    end
    
    fid = fopen(cname, 'r');
    if fid ~= -1
        id = id + 1;
        Objects(id).PropertyNames{1} = 'Name';
        Objects(id).PropertyValues{1} = name;
        Objects(id).PropertyNames{2} = 'SampleRate';
        Objects(id).PropertyValues{2} = double(rate);
        Objects(id).PropertyNames{3} = 'FileDataType';
        Objects(id).PropertyValues{3} = 'int16';
        Objects(id).PropertyNames{4} = 'Type';
        Objects(id).PropertyValues{4} = 'analog';
        Objects(id).PropertyNames{5} = 'Virtual';
        Objects(id).PropertyValues{5} = 0;
        Objects(id).PropertyNames{6} = 'Digitize';
        Objects(id).PropertyValues{6} = 0;
        Objects(id).PropertyNames{7} = 'CustomConversionFactor';     
        Objects(id).PropertyValues{7} = 1;
        if (ifile == 40)
            Objects(id).PropertyValues{7} = 1;   % Accel temperature is in Celcius.
        elseif (ifile >= 36)
            Objects(id).PropertyValues{7} = 0.001;  % Convert Accel axis to gravity.
        elseif (ifile >= 33)
            Objects(id).PropertyValues{7} = 2.4/32768;  % Convert Aux channels to volts.
        elseif (ifile >= 1)
            Objects(id).PropertyValues{7} = 0.195;  % From Intan spec sheet. Converts A/D units to microvolts.
        end
        Objects(id).PropertyNames{8} = 'CustomDataUnits';
        Objects(id).PropertyValues{8} = 'Binary';
        if (ifile == 40)
            Objects(id).PropertyValues{8} = 'C'; % Celcius
        elseif (ifile >= 36)
            Objects(id).PropertyValues{8} = 'g'; % Gravity
        elseif (ifile >= 33)
            Objects(id).PropertyValues{8} = 'V'; % Aux channels in Volts
        elseif (ifile >= 1)
            Objects(id).PropertyValues{8} = 'uV';% Microvolts for biophysical channels
        end
        Objects(id).PropertyNames{9} = 'SignalName';
        Objects(id).PropertyValues{9} = name;
        Objects(id).Parent = 2;
        Objects(id).Children = [];
        Objects(id).Class = 'analog channel';
        Objects(id).SubClass = 'analog signal';
        Objects(id).Changed = 0;
        Objects(id).Version = 2;
        Objects(id).ID = id+1;
        
        filepath = [mdadir name '.chn'];
        outfid = fopen(filepath, 'w');
        if (outfid < 0)
            disp( ['Cound not create: ' filepath]); % Could not open file.
        else
            % Write zero header and binary data
            experimentTime = 0;
            count = fwrite(outfid, zeroHeader, 'integer*2');
            if count == 1024
                % Copy NC3 data file to MDA channel file
                nsamples = 10 * rate; % Buffer size in samples
                if (ifile >= 33) && (ifile <= 35)
                    % Aux channels are unsigned but we right shift by 1 bit
                    % to fit them into int16 values.
                    udata = fread(fid, nsamples, '*uint16');
                    while length(udata) == nsamples
                        fwrite(outfid, bitshift(udata, -1), 'int16');
                        experimentTime = experimentTime + 10;
                        udata = fread(fid, nsamples, '*uint16');
                    end
                    fwrite(outfid, bitshift(udata, -1), 'int16');
                    experimentTime = experimentTime + floor(length(udata) / rate); % Experiment time in whole seconds
                else
                    % Other channels are already int16 compatible.
                    data = fread(fid, nsamples, '*int16');       
                    while length(data) == nsamples
                        fwrite(outfid, data, 'int16');
                        experimentTime = experimentTime + 10;
                        data = fread(fid, nsamples, '*int16');
                    end
                    fwrite(outfid, data, 'int16');
                    experimentTime = experimentTime + floor(length(data) / rate); % Experiment time in whole seconds
                end
                maxTime = max(maxTime, experimentTime);
                disp([experimentName ' Channel ' name ' ' num2str(experimentTime) ' seconds']);
                
                % Write binary header
                fseek(outfid, 0, -1);
                fileID = 16000;
                fwrite(outfid, fileID, 'integer*2');
                dataMark = 2048;
                fwrite(outfid, dataMark, 'integer*2');
                fileVersion = 2;
                fwrite(outfid, fileVersion, 'double');
                fwrite(outfid, experimentTime, 'double');

                % Write text header
                header = sprintf(['SignalName = ''%s''\r\n' ...
                    'Comments = ''''\r\n' ...
                    'UseHardware = 0\r\n' ...
                    'Digitize = 0\r\n' ...
                    'Experiment = ''%s''\r\n' ...
                    'FileDataType = ''int32''\r\n' ...
                    'Name = ''%s''\r\n' ...
                    'SampleRate = %d\r\n' ...
                    'Type = ''analog''\r\n' ...
                    'Virtual = 0\r\n' ...
                    'Changed = 1\r\n' ...
                    'Class = ''analog channel''\r\n' ...
                    'ID = %d\r\n' ...
                    'Subclass = ''analog signal''\r\n' ...
                    'Version = 2\r\n'], name, experimentName, name, rate, id+1);
                fwrite(outfid, header, 'integer*1');               
            else
                disp(['Cound not write to: ' filepath]); % Could not write to file.
            end     
            
            fclose(outfid);  % Close output file.
        end
        fclose(fid); % Close input file
    end % if fid
end % for ifile

% Create Timestamp channels

cname = [fpath fprefix '_Events.u32'];
rate = 5000;
name = 'Events';
fid = fopen(cname, 'r');
if fid ~= -1
    events = fread(fid, [3,inf], 'uint32'); % Read in events.
    fclose(fid);

    for ichan = 0:15

        % Setup for different types of events.
        if (ichan <= 7)
            % Events from discriminators
            name = ['Event' num2str(ichan)];
            index = find(events(1,:) == 1 & (events(2,:) == ichan));
        elseif (ichan <= 13)
            % Stimulation events
            name = 'Stim0';
            istim = ichan - 7; % Stim channel index
            name(end) = name(end) + istim; % StimA..StimC
            index = find(events(1,:) == 2 & (events(2,:) == istim));
        elseif (ichan == 14)
            % Condition changed events
            name = 'Condition';
            index = find(events(1,:) == 0);
        elseif (ichan == 15)
            % Analog sampling status changed events
            name = 'Status';
            index = find(events(1,:) == 3);
        else
            % No event
            index = [];
        end

        if ~isempty(index)
            % Create new timestamp object
            id = id + 1;
            Objects(id).PropertyNames{1} = 'Name';
            Objects(id).PropertyValues{1} = name;
            Objects(id).PropertyNames{2} = 'SampleRate';
            Objects(id).PropertyValues{2} = rate;
            Objects(id).PropertyNames{3} = 'FileDataType';
            Objects(id).PropertyValues{3} = 'uint32';
            Objects(id).PropertyNames{4} = 'Type';
            Objects(id).PropertyValues{4} = 'digital';
            Objects(id).PropertyNames{5} = 'Virtual';
            Objects(id).PropertyValues{5} = 0;
            Objects(id).PropertyNames{6} = 'Digitize';
            Objects(id).PropertyValues{6} = 0;
            if ichan >= 14
                % Condition and Status events have accessory data.
                Objects(id).PropertyNames{7} = 'Events';
                Objects(id).PropertyValues{7} = 1;
                Objects(id).PropertyNames{8} = 'EventFileDataType';
                Objects(id).PropertyValues{8} = 'int8';
            end
            Objects(id).Parent = 2;
            Objects(id).Children = [];
            Objects(id).Class = 'timestamp channel';
            Objects(id).SubClass = 'digital signal';
            Objects(id).Changed = 0;
            Objects(id).Version = 2;
            Objects(id).ID = id+1;

            filepath = [mdadir name '.chn'];
            outfid = fopen(filepath, 'w');
            if outfid < 0
                disp( ['Cound not create: ' filepath]); % Could not open file.
            else
                % Write zero header
                count = fwrite(outfid, zeroHeader, 'integer*2');
                
                % Write timestamps
                timestamps = events(3, index);
                if ichan >= 14
                    % Condition and Status events are tagged with accessory data = Value
                    values = events(2, index);
                    for istamp=1:length(timestamps)
                        fwrite(fid, timestamps(istamp), 'uint32');
                        fwrite(fid, values(istamp), 'char');
                    end
                else
                    % Other channels have no accessory data.
                    fwrite(outfid, timestamps, 'uint32');
                end
                experimentTime = max(maxTime, ceil((timestamps(end)+1) / rate));
                disp([experimentName ' Timestamp ' name ' ' num2str(length(timestamps)) ' events']);

                % Write binary header
                fseek(outfid, 0, -1);
                fileID = 16000;
                fwrite(outfid, fileID, 'integer*2');
                dataMark = 2048;
                fwrite(outfid, dataMark, 'integer*2');
                fileVersion = 2;
                fwrite(outfid, fileVersion, 'double');
                fwrite(outfid, experimentTime, 'double');

                % Write text header
                header = sprintf(['Reference = ''none''\r\n' ...
                    'RefractoryPeriod = 0\r\n' ...
                    'UseHardware = 0\r\n' ...
                    'Comments = ''''\r\n' ...
                    'Digitize = 0\r\n' ...
                    'EventConditionSet = ''(default)''\r\n' ...
                    'Experiment = ''%s''\r\n' ...
                    'FileDataType = ''uint32''\r\n' ...
                    'Filter = 0\r\n' ...
                    'Name = ''%s''\r\n' ...
                    'SampleRate = %d\r\n' ...
                    'SweepFileDataType = ''int16''\r\n' ...
                    'SweepReference = ''''\r\n' ...
                    'SweepSampleRate = 0\r\n' ...
                    'SweepUseHardware = 0\r\n' ...
                    'SweepWindowStart = 0\r\n' ...
                    'SweepWindowStop = 0\r\n' ...
                    'Sweeps = 0\r\n' ...
                    'Type = ''digital''\r\n' ...
                    'Virtual = 0\r\n' ...
                    'Class = ''timestamp channel''\r\n' ...
                    'Changed = 1\r\n' ...
                    'ID = %d\r\n' ...
                    'Subclass = ''digital signal''\r\n' ...
                    'Version = 2\r\n'], experimentName, name, rate, id+1);
                
                if ichan >= 14
                    header = sprintf('%sEventFileDataType = ''int8''\r\nEvents = 1\r\n', header);
                end
                count = fwrite(outfid, header, 'integer*1');
                fseek(outfid, 0, 'eof');
                fclose(outfid);
            end % if outfid
        end % ~isempty(index
    end %for ichan
end % if fid

% Save the experiment file

experimentfname = fullfile(datapath, [experimentName '.mat']);
save(experimentfname, 'Objects', 'Version');

end %function
