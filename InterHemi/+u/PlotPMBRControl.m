clearvars -except SL SLNeuroControl

Cont = struct('C',[],'I',[]); Ipsi = struct('C',[],'I',[]);

% CHANGE PER ANIMAL
% Ubi
Lchn = [31:32]; Rchn = [6,7];

for n = 1:length(SLNeuroControl)
    
    % Everything is left hemi stim for Ubi
    ContRT.C = SLNeuroControl(n).Rbeta(:,Rchn);
    ContRT.I = SLNeuroControl(n).Rbeta(:,Lchn);
    IpsiRT.C = SLNeuroControl(n).Lbeta(:,Rchn);
    IpsiRT.I = SLNeuroControl(n).Lbeta(:,Lchn);
    
    Cpre = struct('C',[],'I',[]); Ipre = struct('C',[],'I',[]);
    Cpre.C = cellfun(@(x) nanmedian(x,2),ContRT.C(1,:),'UniformOutput',0);
    Cpre.I = cellfun(@(x) nanmedian(x,2),ContRT.I(1,:),'UniformOutput',0);
    Ipre.C = cellfun(@(x) nanmedian(x,2),IpsiRT.C(1,:),'UniformOutput',0);
    Ipre.I = cellfun(@(x) nanmedian(x,2),IpsiRT.I(1,:),'UniformOutput',0);
    
    Cont.C{end+1} = cellfun(@(x,y) x-y,ContRT.C(2,:),Cpre.C,'Uniformoutput',0);
    Cont.I{end+1} = cellfun(@(x,y) x-y,ContRT.I(2,:),Cpre.I,'Uniformoutput',0);
    Ipsi.C{end+1} = cellfun(@(x,y) x-y,IpsiRT.C(2,:),Ipre.C,'Uniformoutput',0);
    Ipsi.I{end+1} = cellfun(@(x,y) x-y,IpsiRT.I(2,:),Ipre.I,'Uniformoutput',0);
    
end

Beta = struct('C',[],'I',[]);  %for each hemisphere
Group = struct('C',[],'I',[]);
for i = 1:length(Cont.C)
    temp = Cont.C{i};
    for j = 1:length(temp)
        temp2 = temp{j};
        for k = 1:size(temp2,1)
            Beta.C = [Beta.C,temp2(k,:)];
            Group.C = [Group.C,k*ones(1,size(temp2,2))];
        end
    end
    temp = Ipsi.C{i};
    for j = 1:length(temp)
        temp2 = temp{j};
        for k = 1:size(temp2,1)
            Beta.C = [Beta.C,temp2(k,:)];
            Group.C = [Group.C,(3+k)*ones(1,size(temp2,2))];
        end
    end
    
    temp = Cont.I{i};
    for j = 1:length(temp)
        temp2 = temp{j};
        for k = 1:size(temp2,1)
            Beta.I = [Beta.I,temp2(k,:)];
            Group.I = [Group.I,k*ones(1,size(temp2,2))];
        end
    end
    temp = Ipsi.I{i};
    for j = 1:length(temp)
        temp2 = temp{j};
        for k = 1:size(temp2,1)
            Beta.I = [Beta.I,temp2(k,:)];
            Group.I = [Group.I,(3+k)*ones(1,size(temp2,2))];
        end
    end
    
end

%% Box plots
positions = 1:6;
color1 = [0, 0.4470, 0.7410]; color2 = [0.8500, 0.3250, 0.0980]; % just using default matlab colors
color = [repmat(color1,3,1); repmat(color2,3,1)];
color = flipud(color); 
xt = [1:6];
xtl = {'Pre','Move','Post','Pre','Move','Post'};

figure; subplot(2,1,1);
boxplot(Beta.C,Group.C,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end
title('Control, Beta, Contra Hemi'); 
hold on; xl = xlim;
plot(xl,[0,0],'k--')

c = get(gca, 'Children');
[~,leg] = legend(c([2,5]), 'Contra', 'Ipsi' );
PatchInLegend = findobj(leg, 'type', 'patch');
set(PatchInLegend, 'facealpha', 0.5);
boxplot(Beta.C,Group.C,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
ylim([-3,3]);

subplot(2,1,2);
boxplot(Beta.I,Group.I,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',.5); % set face colors
end
title('Control, Beta, Ipsi Hemi'); 
hold on; xl = xlim;
plot(xl,[0,0],'k--')

boxplot(Beta.I,Group.I,'Notch','on', 'positions', positions,'Whisker',0,'Symbol','');
ylim([-3,3]);





