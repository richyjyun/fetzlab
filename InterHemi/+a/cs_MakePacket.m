%% regenerate SL (takes long time)

%close all, clear all
%u.getMetaData({'Ubi','Igor','Kato'});

%%
close all, clear all, pack
controlfig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 8]);

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'ubi-split_arb415.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

a.MakePacket(SL,fname, controlfig, 1);

%% UBI
close all, clear all, pack
controlfig = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 8]);

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'ubiarb420.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

a.MakePacket(SL,fname, controlfig);

% KATO
%close all, clear all, pack

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'katoarb420.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

a.MakePacket(SL,fname, controlfig);

% IGOR
%close all, clear all, pack

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'igorarb420.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

a.MakePacket(SL,fname, controlfig);
