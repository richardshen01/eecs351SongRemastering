% Matlab code for sharpening song
% Author: Hyugo Weicht

% Takes in a song to be sharpened as well as an impulse response describing
% the sharpening filter:
% output = fullSong*impulseResponse
% Output writes an audiofile w/ name "sharpened.wav"

function a = Sharpening_Filter(song, Fs, impulseResponse)
    % import song & time domain impulse response
    [h,~] = audioread(impulseResponse);

    % zero pad signal
    if (length(h) <= length(song))
        h = [h; zeros(length(song)-length(h),1)];
    
    else
        h = h(1:length(song));
    end
%     stem(abs(fft(h(:,1))))
%     length(h)
%     length(song)
    
    % Perform equation X*H = Y (fourier domain)
    z(:,1) = ifft(fft(song(:,1)).*fft(h));
    z(:,2) = ifft(fft(song(:,2)).*fft(h));
    
    % IFFT + output song
    %audiowrite('sharpened.wav', z, Fs)
    a = z;%convert to a = z
end
