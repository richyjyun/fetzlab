function trig = getKatoTrig(fname,fs,tchn)

train_threshold = 3000;
datfile = [fname '.bin'];
datfid = fopen(datfile, 'r');

if(tchn == 1)
    trig_index = 2;
elseif(tchn == 2)
    trig_index = 19;
end
nchans = 20;

trig_offset = (trig_index-1) * 4;
fseek(datfid, trig_offset, -1);
skip_bytes = (nchans - 1) * 4;
[samples ~] = fread(datfid, '*single', skip_bytes);
fclose(datfid);
trig = find(samples > 0);       % Find triggers
intervals1 = diff([0; trig]);   % Calcualte intervals
trig = (trig(intervals1 > 30) - 1); % Debounce trigger and remove doublets

end