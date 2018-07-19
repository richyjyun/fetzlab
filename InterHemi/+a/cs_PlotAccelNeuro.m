% cs_PlotAccelNeuro

clear; close all; pack;

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
doAccelNeuro(SL);

clear; close all; pack;

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
doAccelNeuro(SL);

clear; close all; pack;

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
doAccelNeuro(SL);

function doAccelNeuro(SL)
fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', ['twitches',char(SL(1).Animal),'.ps']);

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

t = u.PlotAccelNeuro(SL);
for i = 1:length(t)
    print(t(i), '-dpsc2', fname, '-append')
    close(t(i))
end
end


fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', ['AccelNeuroKato.ps']);
% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

for i = 1:length(SL)
    D = str2num(SL(i).Date);
%     if(D>=20170404)
        t = u.PlotAccelNeuro(SL(i));
        print(t, '-dpsc2', fname, '-append')
        close(t)
        
%         if(~strcmp(SL(i).StimHemi,'NaN'))
%         % Plot neural data
%         [left, right, trig, ccep] = u.PlotNeuralERP(SL(i));
%         print(left, '-dpsc2', fname, '-append')
%         close(left)
%         print(right, '-dpsc2', fname, '-append')
%         close(right)
%         print(trig, '-dpsc2', fname, '-append')
%         close(trig)
%         print(ccep, '-dpsc2', fname, '-append')
%         close(ccep)
%     end
        
end







