clear; close all;

date = '20180201';

%% Find correct block
tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blocks = dir(tankpath);
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';
date = date(3:end);
blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));
blockname = blockname(1,:); % only use first 

%% Find times to load
TT = TDT2mat([tankpath,blockname],'TYPE',2);
Dscm = TT.epocs.Dscm;
[val,ind] = findpeaks(Dscm.data); 
ind = ind(val>1000); val = val(val>1000); 
times = [ind(1)-val(1),ind(1)];  
times = Dscm.onset(times);

%% Load
Snips = TDT2mat([tankpath,blockname],'T1',times(1),'T2',times(2),'TYPE',3,'STORE','eNe1'); Snips = Snips.snips.eNe1;

%% Find channels and codes
spkLim = 2000;
chns = []; codes = [];
for chn = 1:96
    for code = 1:3
        ind = Snips.chan == chn & Snips.sortcode == code;
        if(sum(ind) > spkLim)
            chns(end+1) = chn;
            codes(end+1) = code;
        end
    end
end

chnStr = ''; codeStr = '';
for i = 1:length(chns)
    chnStr = sprintf('%s%d/',chnStr,chns(i));
    codeStr = sprintf('%s%d/',codeStr,codes(i));
end
chnStr(end) = []; codeStr(end) = [];

disp(chnStr)
disp(codeStr)