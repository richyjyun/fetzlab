function noise = findNoise(data,fs)

    threshold = 3; % 3 std devs 
    minwindow = 0.05; % window to check for peak
    rmwindow = [-0.05,0.15]; % window to remove around peak
    
    data = abs(zscore(data)); % zscore and rectify for easier handling
    crossings = find(data > threshold); % find threshold crossings
    dcross = diff(crossings); % find first crossings
    rm = find(dcross < minwindow*fs); rm = rm+1;
    crossings(rm) = [];
    
    % find the peaks around threshold crossings
    noise = zeros(1,length(crossings));
    for c = 1:length(crossings)
        if(crossings(c)+round(minwindow*fs) > length(data))
            continue;
        end
        [~,noise(c)] = max(data(crossings(c):crossings(c)+round(minwindow*fs)));
        noise(c) = noise(c) + crossings(c);
    end
    
    % make matrix of indices to remove
    range = round(rmwindow(1)*fs):1:round(rmwindow(2)*fs);
    noise = repmat(noise, length(range), 1) + repmat(range(:), 1, size(noise,2));
    noise(noise<1) = 1; noise(noise>length(data)) = length(data);
        
end