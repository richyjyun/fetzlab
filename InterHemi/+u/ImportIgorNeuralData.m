function [data2, fs, chnm] = ImportIgorNeuralData(fname, fs_target)

if ~exist('fs_target', 'var'), fs_target = 800; end

chL = 'ML'; % name of left hand channel
chR = 'MR'; % name of right hand channel

load([fname '.cfg'],'-mat');
[Nch,chnm] = uichannelsettings(fname);
itarg = strmatch('targets',chnm); if ~UI.plot_fdb, error('no behavior logged to file'); end
irew = strmatch('reward',chnm);
iL = strfind(chnm, chL);
iL = ~cellfun(@isempty, iL, 'unif', 1);
iR = strfind(chnm, chR);
iR = ~cellfun(@isempty, iR, 'unif', 1);
itrig = strmatch('trig',chnm); %if isempty(itrig), error('no trigger channel found'); end

iL = find(iL);
iR = find(iR);

inds = [iL; iR];% these are motor cortex indeces, ordered by hemisphere

chnm = chnm(inds);

fid = fopen([fname '.bin'],'r'); % load data
fseek(fid,0,'eof');
nbytes = ftell(fid); % number of bytes in file
N = nbytes/4/Nch; % number of samples per channel (4 bytes per single-precision sample)
fseek(fid,0,'bof');

downsampleby = round(UI.samprate/fs_target); % take every "downsampleby" samples
fs = UI.samprate/downsampleby;


fseek(fid, (itrig(1)-1)*4, 'bof'); % offset to start of trigger channel
trigdat = fread(fid, N, 'single', 4*(Nch-1));
triggers = find(trigdat/1e3>100); % > 100mV
%triggers = triggers/downsampleby;
clear trigdat

% d = designfilt('bandstopiir','FilterOrder',2, ... % design notch filter
%     'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%     'DesignMethod','butter','SampleRate',fs);

skip = Nch*4*(downsampleby-1);
precision = sprintf('%d*single',Nch);
fseek(fid,0,'bof');
data2 = fread(fid,[1,inf],precision,skip); %[9, inf], 'single')'; % note the transpose to put data in columns
data2 = reshape(data2,[Nch,length(data2)/Nch])';
data2 = data2(:,inds);
fclose(fid);

% for i = 1:length(inds)
%     
%     fseek(fid, (inds(i)-1)*4, 'bof'); % offset to start of lfp channel
%     tmp = fread(fid, N, 'single', 4*(Nch-1));
%     
%     tmp(1:round(UI.samprate*.5)) = tmp(round(UI.samprate*.5)); % delete startup transient
%     
%     %     tmp = OOARB.Utils.RemoveArtifact(tmp, triggers, [-.001 .035], UI.samprate);% remove artifact
%     
%     %     data2(:,i) = decimate(tmp, downsampleby);
%     data2(:,i) = tmp(round(linspace(1,length(tmp),length(tmp)*fs_target/UI.samprate)));
%     
%     %     data2(:,i) = filtfilt(d, data2(:,i));
%     
% end

function [Nch,chnm,chgu] = uichannelsettings(cfgfile)
% function [Nch,chnm,chgu] = uichannelsettings(cfgfile)
%   Read user-interface settings in *.cfg files generated by various
%   gUSBamp programs to determine number and names of channels.
load([cfgfile '.cfg'],'-mat');
Nch = length(find(UI.ch_enabled))+length(find(UI.ga_trigger)); % number of channels
if isfield(UI,'plot_fdb'), Nch = Nch+1; end
if isfield(UI,'plot_fdb') && str2double(cfgfile(end-10:end-3))>=20120202, Nch = Nch+1; end
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
if isfield(UI,'plot_fdb'), chnm{ii+1} = 'targets'; chgu(ii+1) = NaN; end
if isfield(UI,'plot_fdb') && str2double(cfgfile(end-10:end-3))>=20120202, chnm{ii+2} = 'reward'; chgu(ii+2) = NaN; end