% clear
% load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbiFinal.mat');

RP = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');
BETA = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');
GAMMA = figure('position', [0 0 900 1100], 'paperposition', [-.8 -.6 10 12]); set(gcf,'visible','off');

window = [-0.2,1];
n = 1;

for i = 1:length(SLNeuro(n).chID)
    set(0, 'CurrentFigure', RP);
    subplot(8,4,i);
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).LRP(1,i,:)),'linewidth',1.5);
    hold on;
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).RRP(1,i,:)),'linewidth',1.5);
    xlim(window);
    title(SLNeuro(n).chnm{SLNeuro(n).chID(i)});
    
    set(0, 'CurrentFigure', BETA);
    subplot(8,4,i);
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).Lbeta(1,i,:)),'linewidth',1.5);
    hold on;
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).Rbeta(1,i,:)),'linewidth',1.5);
    xlim(window);
    title(SLNeuro(n).chnm{SLNeuro(n).chID(i)});
    
    set(0, 'CurrentFigure', GAMMA);
    subplot(8,4,i);
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).Lgamma(1,i,:)),'linewidth',1.5);
    hold on;
    plot(SLNeuro(n).tneuro,squeeze(SLNeuro(n).Rgamma(1,i,:)),'linewidth',1.5);
    xlim(window);
    title(SLNeuro(n).chnm{SLNeuro(n).chID(i)});
    
end

packet = 'F:\Dropbox\repos\abogaard\efetz\U\PaperFigures\Figure2PacketKato.ps';
print('-painters',RP, '-dpsc2', packet, '-append'); close(RP);
print('-painters',BETA, '-dpsc2', packet, '-append'); close(BETA);
print('-painters',GAMMA, '-dpsc2', packet, '-append'); close(GAMMA);

