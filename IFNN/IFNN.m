clear
clc
close all

%% load in epocs
tankpath = 'S:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blockname = 'Spanky-180605-135243';
nostimtrial = 1;
T1 = 1; T2 = 0;
Epocs = TDT2mat([tankpath,blockname],'T1',T1,'T2',T2,'TYPE',2); Onset = Epocs.epocs.Dscm.onset; Epocs = Epocs.epocs.Dscm.data;
LFPs = TDT2mat([tankpath,blockname],'T1',1,'T2',1.01,'TYPE',4,'STORE','LFPs'); fs = LFPs.streams.LFPs.fs;

%% Search for intervals
dEpocs = cat(1,diff(Epocs),0);
Interends = find(dEpocs < -1000);
for i = 1:length(Interends)
    for j = (Interends(i)-1):-1:1
        if dEpocs(j) < 1
            Interstart(i,1) = j + 1;
            break
        end
    end
end

%% Loading in appropriate variables
vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

% Arrange into a struct. First row is the fields
SL = struct([]);
for i = 2:size(vals,1)
    if(~strcmp(vals(i,1),'Spike Trig'))
        continue;
    end

    cur = length(SL)+1;
    for j = 1:size(vals,2)
        SL(cur).(char(vals(1,j))) = char(vals{i,j});
    end
end

% Find correct block
day = strcat('20',blockname(8:13));

ind = find(strcmp(day,extractfield(SL,'Date')));

% Get channels and corresponding sort codes
Channels = split(SL(ind).Channels,'/');
Codes = split(SL(ind).SortCodes,'/');

chns = []; codes = [];
for i = 1:length(Channels)
    c = split(Codes(i),',');
    for j = 1:length(c)
        chns(end+1) = str2double(Channels(i));
        codes(end+1) = str2double(c(j));
    end
end

%% Load spikes
window = 0.05;
range = round(-window*fs:1:window*fs);

onset = Onset(Interstart(nostimtrial):Interends(nostimtrial));
snips = TDT2mat([tankpath,blockname],'T1',onset(1),'T2',onset(end),'TYPE',3); 
spikes = snips.snips.eNe1.chan;
sortcodes = snips.snips.eNe1.sortcode;
spiketimes = snips.snips.eNe1.ts;

ch = [81,83];
codes = [1,1];

timesteps = round(fs*(onset(end)-onset(1)));
input_events = zeros(3,timesteps);

for i = 1:length(ch)
    findspikes = find(spikes==chns(i) & sortcodes==codes(i));  
    indspiketimes = round(fs*(spiketimes(findspikes)-spiketimes(1)));
    temp = zeros(1,timesteps); temp(indspiketimes) = 1;
    input_events(i,:) = temp;
end

% I load in the whole file, but often only want to run part of the file, so
% I overwrite timesteps.
maxtime = 10; %minutes
timesteps = round(maxtime*60*fs);
input_events = input_events(:,1:timesteps);
clc
%% Network Parameters
nin = length(ch);
nhid = 0;
nout = 1;
nunits = nin + nhid + nout;

ninputs = nin; % Number of input units receiving imposed output spikes.
nhidden = nhid;
psp_slow_decay = 31/32; % Slow decay constant for PSPs. Example = 31/32; 
psp_fast_decay = 7/8;   % fast decay constant for PSPs. Example = 7/8;
depolarization = 400; % Unit hyperpolarization when spike occurs. Example = 10000
fast_potentials = zeros(nunits,1); % Fast decaying part of unit potentials.
slow_potentials = zeros(nunits,1); % Slow decaying part of unit potentials.
ethresh = 400;
thresholds = ethresh * ones(nunits, 1); % Unit potential threshold.  Example = 5000; 
% x is a variable that is required in order to make it so that the weight
% you input for the psp (desiredweight) is the maximum size of the psp.
x = (-log(log(psp_fast_decay)/log(psp_slow_decay)))/(log(psp_fast_decay)-log(psp_slow_decay));
t = 5; % ms (coincident spike interval)
ex = fs*(t/1000); % number of iterations corresponding to 5 ms
frac = psp_slow_decay^ex-psp_fast_decay^ex; % psp amplitude fraction at 5 ms
desiredweight = ethresh/(frac+1); % calculate weight to detect coincident activation
desiredweight = desiredweight/(psp_slow_decay^x-psp_fast_decay^x); % scaling weights so that psp amplitude equals weight
weights = desiredweight * ones(nunits, nunits); % W(j,i) = connection strength from unit(j) to unit(i).
weights(:, 1:nin) = 0; % Input units don't have connections to them (only from them).
weights(end-nout+1:end,:) = 0; % Output units don't have connections to other units.

output_events = zeros(nunits, timesteps);   % Storage for output events.
unit_potentials = zeros(nunits, timesteps); 
for tstep = 1:timesteps
    input_events_v = input_events(:,tstep);

    fast_potentials = psp_fast_decay .* fast_potentials;  % All unit potentials decay
    slow_potentials = psp_slow_decay .* slow_potentials;  % All unit potentials decay

    % Check for PSP delivery from all units -> hidden and output units.
    from_units = double([input_events_v(1:ninputs);output_events(ninputs+1:ninputs+nhidden,tstep);zeros(nout,1)]);
    to_index = (ninputs+1:nunits)'; % Indexes for the non-input units

    incoming_psps =  (from_units' * weights(:, to_index))';   % Matrix multiply [1 x N] * [N x M]
    slow_potentials(to_index) = slow_potentials(to_index) + incoming_psps;
    fast_potentials(to_index) = fast_potentials(to_index) + incoming_psps;

    % Check for units firing when threshold is reached.
    output_events(:,tstep) = (slow_potentials - fast_potentials >= thresholds); % Find units that fire.
    slow_potentials(logical(output_events(:,tstep))) = slow_potentials(logical(output_events(:,tstep))) - depolarization;  % Reset unit potential
    unit_potentials(:, tstep) = slow_potentials - fast_potentials;
    output_events(1:ninputs,tstep) = input_events_v(1:ninputs); % Copy input unit events to output events     
end



