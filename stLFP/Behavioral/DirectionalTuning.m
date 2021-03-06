function [direction,magnitude] = DirectionalTuning(Mani,spk,fs,window,dt,plt)

% 1 is RU (radial (up) is positive) and 2 is FE (flexion (left) is positive)
Mani.data(2,:) = -Mani.data(2,:); % convert for left/right

center = [median(Mani.data(2,:));median(Mani.data(1,:))];

window = 0; % can modify to get an average of the window of manip data rather than the instantaneous
range = -window*fs:window*fs;

spk = round(spk);

direction = nan(length(spk),1);
magnitude = nan(length(spk),1);
% loc = [];
for s = 1:length(spk)
    loc = [mean(Mani.data(2,spk(s)+range));mean(Mani.data(1,spk(s)+range))];
    loc = loc-center;
    direction(s) = (atan2d(loc(2),loc(1)))*pi/180;
    magnitude(s) = sqrt(loc(1)^2+loc(2)^2);
end

end