close all; clear all; pack

vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

% Arrange into a struct. First row is the fields
SL = struct([]);
for i = 2:size(vals,1)
    if(~strcmp(vals(i,1),'Stim Test'))
        continue;
    end
    if(str2num(vals{i,2}) < 20180122)
        continue;
    end
   
    cur = length(SL)+1;
    for j = 1:size(vals,2)
        SL(cur).(char(vals(1,j))) = char(vals{i,j});
    end
end

tankpath = 'Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223\';
blocks = dir('Y:\~NeuroWest\Spanky\SpikeTrigger-180122-105223');
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';
for i = 1:length(SL)
    date = SL(i).Date; date = date(3:end); disp(date);
    blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));
    trigChns = 32;
    stimChn = 32;
    
    for j = 1:size(blockname,1)
        TT = TDT2mat([tankpath,blockname(j,:)],'TYPE',2);
        Dscm = TT.epocs.Dscm;
        [val,ind] = findpeaks(Dscm.data);
        ind = ind(val>1000); val = val(val>1000);
        tests = 1; times = [ind(tests)-val(tests),ind(tests)];
        times = Dscm.onset(times);
        
        SaveAmpPhase(tankpath,blockname(j,:),trigChns,1,stimChn,times);
    end
end


