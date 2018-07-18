function events = nc3events(file_name)
% events = nc3events(file_name)
% Reads event, stimulus, and status timestamps from a Neurochip3 data set.
%
% file_name is partial or full path name to the neurchip3 .mat settings
%   file in the data folder. if omitted, a file dialog box will be
%   displayed to select the .mat file.
%
% events.nevents returns total number of discrim and stim events.  Returns empty on error.
% events.ncondition(1..8) returns number of start time events for each condition.
% events.condition{1..8} returns condition start event times in seconds.
% events.condition_duration(1..8) returns durations for each condition in seconds.
% events.ndiscrim(1..8) returns number of discimination events for each event discriminator.
% events.discrim{1..8} returns discrimination event times in seconds.
% events.nstim(1..6) returns number of stimulations on each stim channel.
% events.stim{1..3} returns stimulation event times in seconds.
% events.missing_sections returns a 2xN matrix of time intervals where
%   analog data doest exist in the analog channel files. Example:
%   events.missing_sections(1,1) is the first unusable time (in seconds).
%   events.missing_sections(2,1) is the next usuable time (in seconds). 
%   if empty, then there are no missing sections of data.
% events.last_event_time returns the time in seconds of the last event
% events.session_time returns duration of the session in seconds.
% events.channel_names{} = returns names for analog channels.
% events.channel_rate() returns sample rates for analog channels.

if nargin < 1
    fname = '';
else
    fname = file_name;
end

% Default return value

events.nevents = [];
events.ncondition = zeros(1,8);
events.condition = {[], [], [], [], [], [], [], []};
events.condition_duration = zeros(1,8);
events.ndiscrim = zeros(1,8);
events.discrim =  {[], [], [], [], [], [], [], []};
events.nstim = zeros(1,6);
events.stim = {[], [], [], [], [], []};
events.nmissing_sections = 0;
events.missing_sections = zeros(2,0);
events.last_event_time = 0;
events.session_time = 0;
events.channel_rate = zeros(1,48);
events.channel_names = cell(1,48);

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

ioffset = 0;
if p.version < 320
    ioffset = 1;  % Adjustment for the zero based condition/discrim/stim indexes in earlier versions.
    events.condition_duration(1:4) = p.stim_cond_duration(1:4); % 4 conditions in v310 and below.
    remote_emg_channels = 0;
else
    events.condition_duration(1:8) = p.stim_cond_duration(1:8); % Up to 8 conditions in v320 and above.
    remote_emg_channels = p.remote_channels;
end

% Grab some information about the channels
events.channel_rate = p.channel_rate; % Channel Sample Rates
for ichan = 1:length(events.channel_rate)
    if (ichan >= 41 + remote_emg_channels)
        events.channel_rate(ichan) = 0; % Cull out unused remote EMG channels
    end
    if events.channel_rate(ichan) > 0
        events.channel_names{ichan} = p.channel_names{ichan}; % Copy Channel Names
        % Grab total session time from first sampled channel or Aux channel
        if (events.session_time == 0) && (ichan <= 35)
            if (ichan <= 32)
                cname = [fpath fprefix '_Chan' num2str(ichan-ioffset,'%.02d') '.i16'];
            else
                cname = [fpath fprefix '_Chan' num2str(ichan-ioffset,'%.02d') '.u16'];
            end
            fid = fopen(cname, 'r');
            if (fid ~= -1)
                fseek(fid, 0, 1); % End of file
                events.session_time = ftell(fid) / (2 * events.channel_rate(ichan)); % Session time in seconds
                fclose(fid);
            end
        end
    end
end

disp(['File version ' num2str(p.version)]);

% Open the events file, and read all events.

cname = [fpath fprefix '_Events.u32'];
fid = fopen(cname, 'r');
if fid == -1
	disp(['error: Events file not found: ' cname]);
    return;
end

all_events = fread(fid, [3, inf], '*uint32')';
fclose(fid);

% Get time of last event in seconds.

events.nevents = 0;
events.last_event_time = double(all_events(end, 3)) / 5000;
if (events.session_time < events.last_event_time)
    events.session_time = events.last_event_time;
end

% Handle conditions

for icond = 1:9
    index = find((all_events(:,1) == 0) & (all_events(:,2) == icond-ioffset));
    events.condition{icond} = double(all_events(index,3)) ./ 5000;
    events.ncondition(icond) = length(index);
end

% Handle events

for idiscrim = 1:8
    index = find((all_events(:,1) == 1) & (all_events(:,2) == idiscrim-ioffset));
    events.discrim{idiscrim} = double(all_events(index,3)) ./ 5000;
    events.nevents = events.nevents + length(index);
    events.ndiscrim(idiscrim) = length(index);
end

% Handle stims

for istim = 1:6
    index = find((all_events(:,1) == 2) & (all_events(:,2) == istim-ioffset));
    events.stim{istim} = double(all_events(index,3)) ./ 5000;
    events.nevents = events.nevents + length(index);
    events.nstim(istim) = length(index);
end

% Handle missing sections

index = find(all_events(:,1) == 3); % find all section markers
n = length(index);
i = 1;
while i <= n
    if all_events(index(i),2) == 0
        % Start of a missing data section
        ts = double(all_events(index(i),3)) / 5000;
        events.missing_sections(:,end+1) = [ts;ts+1];
        i = i + 1;
        while i <= n
           if all_events(index(i),2) == 1
               % end of missing data section
               ts = double(all_events(index(i),3)) / 5000;
               events.missing_sections(2,end) = ts;
               break;
           end
           i = i + 1;
        end
    end
    i = i + 1;
end

% Display file summary info.

events.nmissing_sections = size(events.missing_sections, 2);
missing_time = 0;
if events.nmissing_sections > 0
    missing_time = sum(events.missing_sections(2,:) - events.missing_sections(1,:));
end
disp([num2str(events.nevents) ' events. ' num2str(events.nmissing_sections) ' missing sections of data (' num2str(missing_time) ' seconds combined)']);
