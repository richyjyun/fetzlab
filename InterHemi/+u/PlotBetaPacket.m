% clear;
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorNeuro.mat');
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgorFinal.mat');

packet = 'F:\Dropbox\repos\abogaard\efetz\U\PaperFigures\BetaPacketIgorSmoothed.ps';
window = [-0.5,1.5];
for n = 1:length(SLNeuro)
    
    smth = 50; 
    
    LBETA = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');
    RBETA = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');
    
    for i = 1:length(SLNeuro(n).chID)
        set(0, 'CurrentFigure', LBETA);
        
        subplot(8,4,i);
        sig = squeeze(SLNeuro(n).Lbeta(1,i,:));
        sig = smooth(sig,smth);
        plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
        hold on;
        sig = squeeze(SLNeuro(n).Lbeta(2,i,:));
        sig = smooth(sig,smth);
        plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
        xlim(window);
        title(SLNeuro(n).chnm{SLNeuro(n).chID(i)});
        
        set(0, 'CurrentFigure', RBETA);
        
        subplot(8,4,i);        
        sig = squeeze(SLNeuro(n).Rbeta(1,i,:));
        sig = smooth(sig,smth);
        plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
        hold on;
        sig = squeeze(SLNeuro(n).Rbeta(2,i,:));
        sig = smooth(sig,smth);
        plot(SLNeuro(n).tneuro,sig,'linewidth',1.5);
        xlim(window);
        title(SLNeuro(n).chnm{SLNeuro(n).chID(i)});
        
    end
    
    print('-painters',LBETA, '-dpsc2', packet, '-append'); close(LBETA);
    print('-painters',RBETA, '-dpsc2', packet, '-append'); close(RBETA);
end
