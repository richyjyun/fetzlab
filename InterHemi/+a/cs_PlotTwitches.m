%% Plots twitches based around stimulus trigger
clear all; close all; pack;

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'twitchesKato.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
doTwitch(SL,fname);

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
doTwitch(SL,fname);

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
doTwitch(SL,fname);

function doTwitch(SL,fname)

t = u.PlotTwitches(SL);

for i = 1:length(t)
    print(t(i), '-dpsc2', fname, '-append')
    close(t(i))
end

end
