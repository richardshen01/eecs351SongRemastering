
%gaussian filter whole window, 50Hz high pass
function z = Second_Gaussian(y,Fs)
    w = gausswin(length(y)); % gaussian window
    w = [zeros(length(y) ,1); w; zeros(length(y),1)]; %zero pad signal
    
    
    %perform convolution
    a = fft(w);
    b = fft([zeros(length(y),1); y(:,1); zeros(length(y),1)]);
    c = fft([zeros(length(y),1); y(:,2); zeros(length(y),1)]);
    filteredSignal1 = ifft(a.*b);
    filteredSignal2 = ifft(a.*c);

    z(:,1) = highpass(filteredSignal1(1:length(y)),50,Fs);
    z(:,2) = highpass(filteredSignal2(1:length(y)),50,Fs);
    
    % audiowrite('test.wav',z,Fs);
end