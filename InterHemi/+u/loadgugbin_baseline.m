function [data fs chnm itrigch] = loadgugbin_baseline(fnm, biofb)

% Loads a Guger binary file into Matlab
%
% [data fs channel_names trigger_channels_ids] = loadgugbin(filename, biofb)
%
% filename without suffix (e.g. '2012_01_19_003')
% biofb set to 1 if biofeedback signals are recorded

load([fnm '.cfg'],'-mat');

Nch = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
if biofb==1
    if isfield(UI,'plot_fdb'), Nch = Nch+2;
    else
        warning ('No biofeedback signal detected!')
    end
else
    Nch = Nch+1; % Larry's "extra" trigger-type channel
end
chnm = cell(Nch,1); % name of channels
chgu = zeros(Nch,1); % gUSBamp that channel was recorded on
ii = 0;
for iga = 1:length(UI.ga_trigger)
    ind = find(UI.ch_enabled(iga,:));
    chnm(ii+1:ii+length(ind)) = UI.ch_name(iga,ind)';
    chgu(ii+1:ii+length(ind)) = iga*ones(length(ind),1);
    ii = ii + length(ind);
    if UI.ga_trigger(iga), chnm{ii+1} = ['trig ' num2str(iga)]; chgu(ii+1) = iga; ii = ii + 1; end
end
itrigch = strncmp('trig',chnm,4);
fs = UI.samprate;
sdx = round(15*UI.samprate);
fid = fopen([fnm '.bin'],'r');
eof = 0; s = 1; tr = []; data = [];
while ~eof
    [dat,count] = fread(fid,[Nch,sdx],'single');
    if count < Nch*sdx,
        data(s:s+(count / Nch)-1,:) = dat';
        eof = 1;
        continue;
    else
        data(s:s+sdx-1,:) = dat';
        s = s + sdx;
    end
end
fclose(fid);
end