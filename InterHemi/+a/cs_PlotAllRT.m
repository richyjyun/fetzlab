close all, clear all, pack

fname = fullfile('F:\Dropbox\repos\abogaard\efetz\U\packets', 'RT.ps');

% Delete previous ps file so it doesn't keep appending
if(exist(fname))
    delete(fname);
end

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataUbi.mat')
SL = a.RemoveTwitch(SL);
SL = a.AppendReactionTimes(SL);
SL = a.AppendNormalizedDelay(SL);
[c,p,r] = u.PlotRT(SL);
print(c, '-dpsc2', fname, '-append')
close(c)
print(p, '-dpsc2', fname, '-append')
close(p)
print(r, '-dpsc2', fname, '-append')
close(r)

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataIgor.mat')
% Try to clean up triggers
for i = 1:length(SL)
   if(strcmp(SL(i).StimHemi, 'L'))
       if(any(strmatch('Contra',SL(i).Condition)))
           hasTrig = zeros(1,length(SL(i).righttrials));
           Trials = SL(i).righttrials;
       elseif(any(strmatch('Ipsi',SL(i).Condition)))
           hasTrig = zeros(1,length(SL(i).lefttrials));
           Trials = SL(i).lefttrials;
       else
           continue;
       end
   elseif (strcmp(SL(i).StimHemi, 'R'))
       if(any(strmatch('Contra',SL(i).Condition)))
           hasTrig = zeros(1,length(SL(i).lefttrials));
           Trials = SL(i).lefttrials;
       elseif(any(strmatch('Ipsi',SL(i).Condition)))
           hasTrig = zeros(1,length(SL(i).righttrials));
           Trials = SL(i).righttrials;
       else
           continue;
       end
   else
       continue;
   end
       
   bad = [];
   inds = [];
   
   for j = 1:length(SL(i).trig1)
%        norm = abs(SL(i).trig1(j)-Trials(:,1));
%        ind = find(norm == min(norm),1);
       ind = find(SL(i).trig1(j) > Trials(:,1)-50 & SL(i).trig1(j)<Trials(:,2)+50,1) ;
       if(isempty(ind))
           inds(j) = 0;
           bad(end+1) = j;
           continue;
       end
       inds(j) = ind;
       if(hasTrig(ind))
           bad(end+1) = j;
       else
           hasTrig(ind) = 1;
       end
   end
   
   SL(i).trig1(bad) = [];
   
end
SL = a.RemoveTwitch(SL);
SL = a.AppendReactionTimes(SL);
SL = a.AppendNormalizedDelay(SL);
[c,p,r] = u.PlotRT(SL);
print(c, '-dpsc2', fname, '-append')
close(c)
print(p, '-dpsc2', fname, '-append')
close(p)
print(r, '-dpsc2', fname, '-append')
close(r)

load('F:\Dropbox\repos\abogaard\efetz\U\Code\+u\MetaDataKato.mat')
SL = a.RemoveTwitch(SL);
SL = a.AppendReactionTimes(SL);
SL = a.AppendNormalizedDelay(SL);
[c,p,r] = u.PlotRT(SL);
print(c, '-dpsc2', fname, '-append')
close(c)
print(p, '-dpsc2', fname, '-append')
close(p)
print(r, '-dpsc2', fname, '-append')
close(r)

