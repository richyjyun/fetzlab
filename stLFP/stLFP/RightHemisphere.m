close all; clear all; pack

vals = GetGoogleSpreadsheet('1WLfx_3Zq1MdA2T0S6-LUTU0QqMe68vDLM8caA5_lEMc');

% Arrange into a struct. First row is the fields
SL = struct([]);
for i = 2:size(vals,1)
    if(~strcmp(vals(i,1),'Stim Test'))
        continue;
    end

    cur = length(SL)+1;
    for j = 1:size(vals,2)
        SL(cur).(char(vals(1,j))) = char(vals{i,j});
    end
end

tankpath = 'Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939';
blocks = dir(tankpath);
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';
for i = 1:length(SL)    
    date = SL(i).Date; date = date(3:end); disp(date);
    
    if(str2num(date)~=180116), continue; end;
    
    fname = fullfile('F:\S\Packets', [date,'.ps']);
    % Delete previous ps file so it doesn't keep appending
    if(exist(fname))
        continue;
    end

    blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));
    if(size(blockname,1)>1)
        continue;
    end
    % Get channels and corresponding sort codes
    Channels = split(SL(i).Channels,'/');
    Codes = split(SL(i).SortCodes,'/');
    
    chns = []; codes = [];
    for i = 1:length(Channels)
        c = split(Codes(i),',');
        for j = 1:length(c)
            chns(end+1) = str2double(Channels(i));
            codes(end+1) = str2double(c(j));
        end
    end
    
    TT = TDT2mat([tankpath,'\',blockname],'TYPE',2);
    Dscm = TT.epocs.Dscm;
    [val,ind] = findpeaks(Dscm.data);
    ind = ind(val>5000); val = val(val>5000);
    times = [ind(1)-val(1),ind(1)];
    times = Dscm.onset(times)';
        
    PlotWstLFP(blockname,chns,codes,times,fname,0);
end