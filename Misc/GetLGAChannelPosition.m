function [p,c] = GetLGAChannelPosition(channel)
% Takes in a channel number and returns position on LGA pad (p) in a 19x20
% subplot. Also returns the channel name (c)
%
% RJY 08/27/18

% Set up plot positions
pos1 = 2:2:18; pos1 = [pos1',pos1'+1];
pos2 = 1:2:19; pos2 = [pos2',pos2'+1];
pos = [];
for i = 1:19 % 19 rows in LGA pad
    temp = (i-1)*20; % 20 columns
    align = mod(i,2);
    if(align)
        temp = temp+pos1;
    else
        temp = temp+pos2;
    end
    pos = [pos;temp];
end

% Define positions that will be used
good = [{4:6};
    {4:7};
    {3:7};
    {3:8};
    {2:8};
    {2:9};
    {1:9};
    {1:10};
    {1:9};
    {1:10};
    {1:9};
    {1:10};
    {1:9};
    {2:9};
    {2:8};
    {3:8};
    {3:7};
    {4:7};
    {4:6}];
shift = 0;
for i = 1:length(good)
    good{i} = good{i}+shift;
    align = mod(i,2);
    if(align)
        shift = shift+length(pos1);
    else
        shift = shift+length(pos2);
    end
end
good = cell2mat(good');
pos = pos(good,:);

% Channel names as the order on the pad
chnnm = [{'Ref'},{'Gnd'},{'Ref'},...
        {'C32'},{'C31'},{'D32'},{'D30'},...
        {'C29'},{'B31'},{'A32'},{'D31'},{'D28'},...
        {'C28'},{'C30'},{'B29'},{'A30'},{'D29'},{'D26'},...
        {'C26'},{'C27'},{'B27'},{'B32'},{'A28'},{'D27'},{'D24'},...
        {'C24'},{'C25'},{'B25'},{'B28'},{'A31'},{'A26'},{'D25'},{'D22'},...
        {'C21'},{'C23'},{'B23'},{'B26'},{'B30'},{'A29'},{'A24'},{'D23'},{'D20'},...
        {'C19'},{'C22'},{'B21'},{'B24'},{'B20'},{'A21'},{'A27'},{'A22'},{'D21'},{'D18'},...
        {'C20'},{'B17'},{'B19'},{'B14'},{'B22'},{'A25'},{'A20'},{'A18'},{'D19'},...
        {'C17'},{'C18'},{'B15'},{'B12'},{'B18'},{'A19'},{'A23'},{'A16'},{'D17'},{'D16'},...
        {'C16'},{'B13'},{'B11'},{'B10'},{'A15'},{'A13'},{'A12'},{'A14'},{'D15'},...
        {'C15'},{'C14'},{'B9'},{'B8'},{'B16'},{'A17'},{'A11'},{'A10'},{'D13'},{'D14'},...
        {'C13'},{'C12'},{'B7'},{'B6'},{'A5'},{'A9'},{'A8'},{'D11'},{'D12'},...
        {'C11'},{'C10'},{'B5'},{'B4'},{'A7'},{'A6'},{'D9'},{'D10'},...
        {'C9'},{'C8'},{'B3'},{'B2'},{'A4'},{'D7'},{'D8'},...
        {'C7'},{'C6'},{'B1'},{'A3'},{'D5'},{'D6'},...
        {'C5'},{'C4'},{'A1'},{'D2'},{'D4'},...
        {'C3'},{'C2'},{'A2'},{'D3'},...
        {'C1'},{'Gnd'},{'D1'}];

% convert to numbers
chnnum = zeros(1,length(chnnm));
for i = 1:length(chnnm)
    if(strcmp(chnnm{i},'Ref'))
        chnnum(i) = -1;
    elseif(strcmp(chnnm{i},'Gnd'))
        chnnum(i) = 0;
    else
        shift = 0;
        switch chnnm{i}(1)
            case 'A'
                shift = 0;
            case 'B' 
                shift = 32;
            case 'C'
                shift = 64;
            case 'D'
                shift = 96;
        end
        chnnum(i) = shift + str2num(chnnm{i}(2:end));
    end
end

% Return appropriate position and channel name
ind = find(chnnum==channel);
if(isempty(ind))
    error('Invalid channel number')
else
    p = pos(ind,:);
    c = chnnm{ind};
end

end


