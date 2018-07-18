close all; clear all; pack

fname = fullfile('F:\S\Packets', 'Whitened.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

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

blocks = dir('Y:\~NeuroWest\Spanky\CycleTriggered-170710-143939');
blocks = blocks([blocks.isdir]);
blocks = extractfield(blocks,'name')';
for i = 1:length(SL)
    date = SL(i).Date; date = date(3:end); disp(date);
    blockname = char(blocks(find(~cellfun(@isempty,strfind(blocks,date)))));
    trigChns = str2double(SL(i).Trigger);
    stimChn = str2double(SL(i).Stim);
    times = double(split(SL(i).Tests,'/'));
    times = reshape(times,2,length(times)/2)';
    [allChns,avgDist] = AdapterTest(blockname,trigChns,stimChn,times);
    print(allChns, '-dpsc2', fname, '-append');
    print(avgDist, '-dpsc2', fname, '-append');
    close(allChns); close(avgDist);
end