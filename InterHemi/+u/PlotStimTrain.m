pw = 0.2 ; period = 30;
pulses = 10;
shift = period * 5;

x = [0,pw,pw*2]; 

for i = 1:pulses - 1
    x = [x,x(end-2:end)+period];
end

x = sort([x,x]); x = [0,x+shift,x(end)+shift*2];

y = [0,1,1,-1,-1,0]; y = repmat(y,1,pulses); y = [0,y,0];

figure; plot(x,y,'k'); yticks(''); xticks([shift:60:shift+270]); xticklabels({'0','60','120','180','240'})
xlabel ('ms'); box off; ylim([-1.5,1.5])