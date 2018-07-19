% Run u.PlotBetaGamma or u.PlotBetaRebound first
% For ubi (cond & post, rebound) - Gamma goes up the slower the RT is. Stronger for hemisphere
% during ipsilateral trial. Nothing in beta. Hemisphere is more active,
% limiting the contralateral hemisphere from reacting faster???

figure;
Rchns = 1:16; Lchns = 17:32;

% left hemisphere
% left trials
for i = Lchns
    subplot(2,2,1); hold on;
    scatter(RT.L,Beta.L(:,i),'k.'); title('Left Trial'); ylabel('Beta');
    subplot(2,2,3); hold on;
    scatter(RT.L,Gamma.L(:,i),'k.'); ylabel('Gamma');
end
LagsL = repmat(RT.L,1,size(Beta.L,2)/2);
BetaL = reshape(Beta.L(:,Lchns),1,numel(Beta.L)/2);
GammaL = reshape(Gamma.L(:,Lchns),1,numel(Gamma.L)/2);
tbl = table(LagsL',BetaL',GammaL','VariableNames',{'Lag','Beta','Gamma'});
fitLB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitLG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
LBetaL = BetaL;

%right trials
for i = Lchns
    subplot(2,2,2); hold on;
    scatter(RT.R,Beta.R(:,i),'k.'); title('Right Trial');
    subplot(2,2,4); hold on;
    scatter(RT.R,Gamma.R(:,i),'k.');
end
LagsR = repmat(RT.R,1,size(Beta.R,2)/2);
BetaR = reshape(Beta.R(:,Lchns),1,numel(Beta.R)/2);
GammaR = reshape(Gamma.R(:,Lchns),1,numel(Gamma.R)/2);
tbl = table(LagsR',BetaR',GammaR','VariableNames',{'Lag','Beta','Gamma'});
fitRB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitRG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
LBetaR = BetaR;

figure;
% right hemisphere
% left trials
for i = Rchns
    subplot(2,2,1); hold on;
    scatter(RT.L,Beta.L(:,i),'k.'); title('Left Trial'); ylabel('Beta');
    subplot(2,2,3); hold on;
    scatter(RT.L,Gamma.L(:,i),'k.'); ylabel('Gamma');
end
LagsL = repmat(RT.L,1,size(Beta.L,2)/2);
BetaL = reshape(Beta.L(:,Rchns),1,numel(Beta.L)/2);
GammaL = reshape(Gamma.L(:,Rchns),1,numel(Gamma.L)/2);
tbl = table(LagsL',BetaL',GammaL','VariableNames',{'Lag','Beta','Gamma'});
fitLB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitLG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
RBetaL = BetaL;

% right trials
for i = Rchns
    subplot(2,2,2); hold on;
    scatter(RT.R,Beta.R(:,i),'k.'); title('Right Trial');
    subplot(2,2,4); hold on;
    scatter(RT.R,Gamma.R(:,i),'k.');
end
LagsR = repmat(RT.R,1,size(Beta.R,2)/2);
BetaR = reshape(Beta.R(:,Rchns),1,numel(Beta.R)/2);
GammaR = reshape(Gamma.R(:,Rchns),1,numel(Gamma.R)/2);
tbl = table(LagsR',BetaR',GammaR','VariableNames',{'Lag','Beta','Gamma'});
fitRB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitRG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
RBetaR = BetaR;
% 
% Lerror = std(BetaL)/sqrt(length(BetaL)); Rerror = std(BetaR)/sqrt(length(BetaR));
% remove = (BetaL > mean(BetaL)+Lerror*2 | BetaL < mean(BetaL)-Lerror*2 | BetaR > mean(BetaR)+Rerror*2 | BetaR < mean(BetaR)-Rerror*2);
% remove = hampel(BetaL)
% [bh,bp] = ttest(BetaL(~remove),BetaR(~remove));
% figure; boxplot([BetaL',BetaR'])
% 
% Lerror = std(GammaL)/sqrt(length(GammaL)); Rerror = std(GammaR)/sqrt(length(GammaR));
% remove = (GammaL > mean(GammaL)+Lerror*2 | GammaL < mean(GammaL)-Lerror*2 | GammaR > mean(GammaR)+Rerror*2 | GammaR < mean(GammaR)-Rerror*2);
% [gh,gp] = ttest(GammaL(~remove),GammaR(~remove));
% figure; boxplot([GammaL',GammaR'])

Lchns = 1:(size(Beta.L,2)/2);
Rchns = Lchns + (size(Beta.L,2)/2);
Rchns = 1:16; Lchns = 17:32;
RBetaL = Beta.L(:,Rchns); RBetaL = RBetaL(:)';
RBetaR = Beta.R(:,Rchns); RBetaR = RBetaR(:)';
LBetaR = Beta.R(:,Lchns); LBetaR = LBetaR(:)';
LBetaL = Beta.L(:,Lchns); LBetaL = LBetaL(:)';


[~,j] =hampel(LBetaL);
HLBetaL = LBetaL(~j);
[~,j] =hampel(LBetaR);
HLBetaR = LBetaR(~j);
[~,j] =hampel(RBetaL);
HRBetaL = RBetaL(~j);
% HRLagsL = RLagsL(~j);
[~,j] =hampel(RBetaR);
HRBetaR = RBetaR(~j);

% tbl = table(HRLagsL',HRBetaL','VariableNames',{'Lag','Beta'});
% fitRL = fitlm(tbl,'Beta~Lag','RobustOpts','on');
% 
% AllLags = unique(HRLagsL);
% BetaChange = zeros(1,length(AllLags));
% for i = 1:length(AllLags)
%     BetaChange(i) = median(HRBetaL(HRLagsL==AllLags(i)));
% end
% figure
% plot(AllLags,BetaChange)

AllBeta = [HRBetaR,HRBetaL,HLBetaR,HLBetaL];
groups = [zeros(1,length(HRBetaR)),ones(1,length(HRBetaL)),2*ones(1,length(HLBetaR)),3*ones(1,length(HLBetaL))];
figure; boxplot(AllBeta,groups); xticklabels([]); box off; ylabel('\Delta PMBR')


% AllBeta = [HLBetaL,HLBetaR,HRBetaL,HRBetaR];
% groups = [zeros(1,length(HLBetaL)),ones(1,length(HLBetaR)),2*ones(1,length(HRBetaL)),3*ones(1,length(HRBetaR))];
% figure; boxplot(AllBeta,groups); xticklabels([]); box off; ylabel('\Delta PMBR')


% [p,t,stats] = anova1(AllBeta,groups);
% [c,m,h,nms] = multcompare(stats);



%% Plot gamma for ipsilateral side
% left hemisphere
% left trials
figure; subplot(2,1,1);
h = plot(fitLG); h(1).Marker = '.'; h(1).Color = [0,0,0]; legend HIDE; h(2).LineWidth = 2; h(3).LineWidth = 1.5; h(4).LineWidth = 1.5;
title('Left Hemisphere, Left Trials'); ylabel('\Delta Gamma Power','Interpreter','tex'); xlabel(''); hold on;
h = plot(nan,nan,'w'); legend(h,['Slope: ',num2str(fitLG.Coefficients{2,1})]); legend('boxoff')
ylim([-1,1]); xlim([-30,25])

subplot(2,1,2);
h = plot(fitRG); h(1).Marker = '.'; h(1).Color = [0,0,0]; legend HIDE; h(2).LineWidth = 2; h(3).LineWidth = 1.5; h(4).LineWidth = 1.5;
title('Right Hemisphere, Right Trials'); ylabel('\Delta Gamma Power','Interpreter','tex');
xlabel('\Delta RT','Interpreter','tex'); hold on;
h = plot(nan,nan,'w'); legend(h,['Slope: ',num2str(fitRG.Coefficients{2,1})]); legend('boxoff')
ylim([-1,1]); xlim([-35,30])


%% For Igor (ipsi vs contra)

figure;
Rchns = 10:18; Lchns = 1:9;

% left hemisphere
% ipsi trials
Lags = []; BetaI = []; GammaI = [];
for j = 1:length(StimHemi)
    if(strcmp(StimHemi(j),'L'))
        chns = Lchns;
    else
        chns = Rchns;
    end
    subplot(2,2,1); hold on;
    scatter(repmat(RT.L(j),1,length(chns)),Beta.L(j,chns),'k.'); title('Ipsi Trial'); ylabel('Beta');
    subplot(2,2,3); hold on;
    scatter(repmat(RT.L(j),1,length(chns)),Gamma.L(j,chns),'k.'); ylabel('Gamma');
    Lags = [Lags,repmat(RT.L(j),1,length(chns))]; BetaI = [BetaI,Beta.L(j,chns)]; GammaI = [GammaI,Gamma.L(j,chns)];
end

tbl = table(Lags',BetaI',GammaI','VariableNames',{'Lag','Beta','Gamma'});
fitIB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitIG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
LBetaL = BetaL;

% left hemisphere
% contra trials
Lags = []; BetaC = []; GammaC = [];
for j = 1:length(StimHemi)
    if(strcmp(StimHemi(j),'L'))
        chns = Rchns;
    else
        chns = Lchns;
    end
    subplot(2,2,1); hold on;
    scatter(repmat(RT.L(j),1,length(chns)),Beta.L(j,chns),'k.'); title('Ipsi Trial'); ylabel('Beta');
    subplot(2,2,3); hold on;
    scatter(repmat(RT.L(j),1,length(chns)),Gamma.L(j,chns),'k.'); ylabel('Gamma');
    Lags = [Lags,repmat(RT.L(j),1,length(chns))]; BetaC = [BetaC,Beta.L(j,chns)]; GammaC = [GammaC,Gamma.L(j,chns)];
end

tbl = table(Lags',BetaC',GammaC','VariableNames',{'Lag','Beta','Gamma'});
fitIB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitIG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');


figure;
% right hemisphere
% left trials
for i = Rchns
    subplot(2,2,1); hold on;
    scatter(RT.L,Beta.L(:,i),'k.'); title('Left Trial'); ylabel('Beta');
    subplot(2,2,3); hold on;
    scatter(RT.L,Gamma.L(:,i),'k.'); ylabel('Gamma');
end
LagsL = repmat(RT.L,1,size(Beta.L,2)/2);
BetaL = reshape(Beta.L(:,Rchns),1,numel(Beta.L)/2);
GammaL = reshape(Gamma.L(:,Rchns),1,numel(Gamma.L)/2);
tbl = table(LagsL',BetaL',GammaL','VariableNames',{'Lag','Beta','Gamma'});
fitLB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitLG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
RBetaL = BetaL;

% right trials
for i = Rchns
    subplot(2,2,2); hold on;
    scatter(RT.R,Beta.R(:,i),'k.'); title('Right Trial');
    subplot(2,2,4); hold on;
    scatter(RT.R,Gamma.R(:,i),'k.');
end
LagsR = repmat(RT.R,1,size(Beta.R,2)/2);
BetaR = reshape(Beta.R(:,Rchns),1,numel(Beta.R)/2);
GammaR = reshape(Gamma.R(:,Rchns),1,numel(Gamma.R)/2);
tbl = table(LagsR',BetaR',GammaR','VariableNames',{'Lag','Beta','Gamma'});
fitRB = fitlm(tbl,'Beta~Lag','RobustOpts','on');
fitRG = fitlm(tbl,'Gamma~Lag','RobustOpts','on');
RBetaR = BetaR;