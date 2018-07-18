tankpath = 'F:\S\Test Data\CycleTriggered-170710-143939';
blockname = 'Spanky-170712-143224';

TT = actxcontrol('TTank.X');
TT.ConnectServer('Local','Me');
TT.OpenTank(tankpath,'r');
TT.SelectBlock(blockname);

TT.SetGlobalV('T1', 1000);
TT.SetGlobalV('T2', 1060);

SUA = TT.ReadWavesV('SUA_');

% read first channel
TT.SetGlobalV('Channel', 1);
LFP = TT.ReadWavesV('LFP_');
% preallocate big array
nchan = 96;
LFP = [LFP zeros(length(LFP), nchan-1)];
% read the rest of the channels
for i = 2:nchan
    TT.SetGlobalV('Channel', i);
    LFP(:,i) = TT.ReadWavesV('LFP_');
end

Discrim = TT.ReadEventsSimple('StS1');
Snips = TT.ParseEvV(0, Discrim);
t = TT.ParseEvInfoV(0, Discrim, 6);

TT.CreateEpocIndexing;
Stim = TT.ReadEventsSimple('Stim');
timestamps = TT.ParseEvInfoV(0, Stim, 6);




TT.CloseTank;
TT.ReleaseServer;