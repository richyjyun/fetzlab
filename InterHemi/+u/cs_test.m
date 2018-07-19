clear;

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiAligned.mat')
SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); 
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat','SL','-v7.3');
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); 
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat','SL','-v7.3');
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoAligned.mat')
SL = a.RemoveTwitch(SL); SL = u.AppendThreshold(SL); SL = a.AppendRT(SL); SL = a.AppendNormalizedDelay(SL); 
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKatoFinal.mat','SL','-v7.3');


clear;

dbstop if error
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');
tic;
SLNeuro = u.getLFP(SL);
toc
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiNeuro.mat','SLNeuro','-v7.3');

clear;

dbstop if error
load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');
tic;
SLNeuro = u.getLFP(SL);
toc
save('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorNeuro.mat','SLNeuro','-v7.3');