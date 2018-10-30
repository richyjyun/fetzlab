clear; close all;

%% Load in proper times
% tankpath = 'Y:\~NeuroWest\Spanky\RandomStim-180314-124242\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
% tankpath = 'Y:\~NeuroWest\Spanky\IFNN\';
tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939\';
% tankpath = 'Y:\~NeuroWest\Spanky\Connectivity-180207-131758\';
% tankpath = 'Y:\~NeuroWest\Spanky\SpankyUpstairs\';
blockname = 'Spanky-170921-141010';
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
val(end+1) = Dscm.data(end); ind(end+1) = length(Dscm.data);
ind = ind(val>1000); val = val(val>1000); 
tests = [1,3];
times = [ind(tests)-val(tests),ind(tests)];%-val(tests)+val(2)];  % times = [ind(1) - val(1),ind(1)-val(1)+15000;ind(3)-val(3),ind(3)-val(3)+10000];
times(times==0) = 1;
times = Dscm.onset(times);
% 
% times = [ind(1)-10000,ind(1);ind(3)-val(3),ind(3)-val(3)+10000];

%% Get channels and corresponding sort codes
vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

date = ['20',blockname(8:13)];

% Arrange into a struct. First row is the fields
SL = struct;
for i = 2:size(vals,1)
    if(strcmp(vals(i,2),date))
        for j = 1:size(vals,2)
            SL.(char(vals(1,j))) = char(vals{i,j});
        end
    end
end

Channels = split(SL.Channels,'/');
Codes = split(SL.SortCodes,'/');

chns = []; codes = [];
for i = 1:length(Channels)
    c = split(Codes(i),',');
    for j = 1:length(c)
        chns(end+1) = str2double(Channels(i));
        codes(end+1) = str2double(c(j));
    end
end

%% For cleaned stlfp code
% chns = '11/28/30/32/45/47/53/69/71/77/79/81/82/83/84/88/89/90/93/96';
% chns = strsplit(chns,'/');
% chns = cellfun(@str2num,chns); 
% chns = [45,90,93];
% blockname = ['Spanky-180815-120328';'Spanky-180815-125747'];
tic;
stLFP('tankpath',tankpath,'blockname',blockname,'trigChns',chns,'codes',codes,...
    'times',times,'stimChn',83,'filt',[15,25],'packet',true,'window',[-0.05,0.05],'z',true)
toc;


%% For stlfp
date = blockname(8:13);
fname = fullfile('F:\S\Packets', [date,'.ps']);
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end
% fig = figure;%('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); %set(gcf,'visible','off');;
tic;
PlotWstLFP(tankpath,blockname,30,1,times,fname,0)
toc;

%% For stlfp over time
date = blockname(8:13);
fname = fullfile('F:\S\Packets', [date,'Time.ps']);
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end
tic;
PlotOverTime(tankpath, blockname,32,1,times,fname,[32,60,64])
toc;

%% For Spectra
date = blockname(8:13);
fname = fullfile('F:\S\Packets', [date,'Spectra.ps']);
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end
tic;
SpectralAnalysis(tankpath,blockname,32,1,64,times,fname)
toc;

%% For phase
date = blockname(8:13);
fname = fullfile('F:\S\Packets', [date,'Phase.ps']);
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end
tic;
plotPhaseHist(tankpath,blockname,32,1,times,fname)
toc;
