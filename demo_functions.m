% main code for generating gaussian filtered samples
% current functions:
%   - Main body function
%   - Gaussian filter
%   - Moving Average filter
%   - low pass filter (non standard definiton, read function header)
%   - weakest frequency removal filter (read function header)


% NOTE: MAKE SURE N%4 = 0 FOR INDEXING REASONS

% >takes in the name of a song and outputs a filtered song of length len in
% seconds along with a sample of the original with length len for a direct 
% comparison. 
% >nameIn is input song name
% >nameOut is output name of samples to save. Will save two samples: The 
% >shortened original of form "nameOut_original.wav" , and the filtered 
% sample of form "nameOut_type.wav", where type is gaussian, movAv, etc
% l is length of song sample(in seconds)
% >windowSize is size of filtering window
% >mode is a number between 1 and 2 that has the following options:
%   1 - moving average filter
%   2 - gaussian filter
% >additionalFilter is number between 1 and 3 that has the following options:
%   1 - no additional filtering beyond chosen selection
%   2 - low pass filter to keep percent chosen in percentCutoff
%   3 - remove <percentCutoff>% of weakest amplitude frequencies from 
%       frequency domain 
% >percentCutoff is number between 0 and 1 that corresponds to percent kept
% in low pass filter or percent of weakest frequencies removed for the
% weakest amplitude removal filter
% >LP_dB_Reduction refers to the attenuation performed by low pass filters, 
% where dbReduction = -20 -> 0.1 gain,
% >HPCO is the high pass cutoff frequency in Hz. Recomended value is
% somewhere in ballpark of 100 Hz +/- 50 Hz to remove DC offset

function [a, Fs] = demo_functions(nameIn, mode, windowSize, additionalFiltering)
    [y,Fs] = audioread(nameIn);
    samples = [1, 30*Fs];
    clear y Fs;
    [y,Fs] = audioread(nameIn,samples);
    
    %moving average
    if (mode == 1)
        disp("starting moving average filtering")
        z = movAv(y, windowSize);
        %Homomorphic filtering
        z = decon_function(z, Fs, 1);
        %additional filtering
        if(additionalFiltering == 2)
            z = lowPass(z, 0.8, -20);
        elseif (additionalFiltering == 3)
            z = lowestRemoval(z, 0.8, -20);
        end
        
        %output name and file
%         temp = join([string(nameOut),"_windowSize",string(windowSize), "_movAv.wav"]);
%         audiowrite(temp, z, Fs);
    
    %gaussian
    elseif (mode == 2)
        disp("starting gaussian filtering")
        z = gaussian(y, windowSize, Fs, 100);
        %Homomorphic filtering
        z = decon_function(z, Fs, 1);
        %additional filtering
        if(additionalFiltering == 2)
            z = lowPass(z, 0.8, -20);
        elseif (additionalFiltering == 3)
            z = lowestRemoval(z, 0.8, -20);
        end
        
        %output name and file
%         temp = join([string(nameOut),"_windowSize",string(windowSize),"_gaussian.wav"]);
%         audiowrite(temp, z, Fs);
    end
    a = z; %change to z
    disp("done")
    disp('length of output file')
    disp(length(z))
end

% takes in an input sequence y and outputs a causal moving average of the
% sequence z, with window size windowSize
function z = movAv(y, windowSize)
    w = (1/windowSize)*ones(1,windowSize);
    
    filteredSignal1 = conv(y(:,1),w);
    filteredSignal2 = conv(y(:,2),w);

    z(:,1) = filteredSignal1;
    z(:,2) = filteredSignal2;
end


% takes in an input sequence y and outputs a gaussian filtered sequence z
% assumes that N < length(y) & N is divisible by 4
function z = gaussian(y, N, Fs, CO)
    w = gausswin(N); % gaussian window
    w = [zeros(length(y) ,1); w; zeros(length(y)-N,1); zeros(length(y),1)]; %zero pad signal
    z = double.empty(0,2);
    
    %shift right by 1/4N
    y = circshift(y,N/4);
    
    %loop thru and sum unmixed portions of song
    for i = 1:(((((N + length(y) - mod(length(y),N))/N)))*2)
        %fft
        a = fft(w);
        b = fft([zeros(length(y),1); y(:,1); zeros(length(y),1)]);
        c = fft([zeros(length(y),1); y(:,2); zeros(length(y),1)]);
        filteredSignal1 = ifft(a.*b);
        filteredSignal2 = ifft(a.*c);
        
        %temporary values
        temp1(:,1) = filteredSignal1;
        temp1(:,2) = filteredSignal2;
        temp2 = z;
        
        preHP = temp1(((length(temp1)-length(y)+1))+ N/4: ((length(temp1) - length(y)+N))- N/4,:);
        
        %append the relevant portion of the filtered signal to output sig.
        z = vertcat(temp2 , highpass(preHP,CO,Fs));
        %shift y around by the length of the gaussian window    
        
        y = circshift(y,-N/2);
        %z = vertcat(temp1(1:(N),:),temp2);
    end
    %trim off any excess length
    z = z(1:(length(z) - N + mod(length(y),N)),:);
end

% This function acts as a low pass filter for a given proportion of the
% signal 
% y is the input signal, cutoffPercent is percent of total frequency to
% keep, i.e. if you want to keep lower half of all frequency content,
% cutoffPercent = 0.5
% dB is the attenuation gain of frequencies in the attenuation band
function z = lowPass(y, cutoffPercent, dB)
    %plot, remove comments to enable
    %subplot(2,1,1)
    %stem(abs(fft(y(:,1))))
    %set(gca,'yscal','log')
    %title('input')
    
    %dB to gain
    gain = 10^(dB/20);
    
    %take fft of input signal
    Y1 = fft(y(:,1));
    Y2 = fft(y(:,2));
    
    %determine the indexes corresponding to cutoffPercent
    cutoffIndex1 = uint64(length(Y1)*(cutoffPercent*.5));
    cutoffIndex2 = uint64(length(Y1)*(1-(cutoffPercent*.5)));
    
    for k = cutoffIndex1:cutoffIndex2
        Y1(k) = gain*Y1(k);
        Y2(k) = gain*Y2(k);
    end
    
    z(:,1) = ifft(Y1);
    z(:,2) = ifft(Y2);
   
    %plot, remove comments to enable
    %subplot(2,1,2)
    %stem(abs(fft(z(:,1))))
    %set(gca,'yscal','log')
    %title('output')
end

% This function eliminates X% of the weakest frequencies present in the 
% Fourier domain
% cutoffPercent is the percent of frequencies removed, so for example if
% you wanted to remove the weakest 15% of frequencies, cutoffPercent = 0.15
% dB is the attenuation gain of frequencies in the attenuation band
function z = lowestRemoval(y, cutoffPercent, dB)
    %plot, remove comments to enable
    %subplot(2,1,1)
    %stem(abs(fft(y(:,1))))
    %set(gca,'yscal','log')
    %title('input')

    %dB to gain
    gain = 10^(dB/20);
    
    %total FFT
    Y1 = fft(y(:,1));
    Y2 = fft(y(:,2));

    %magnitude spectrum
    ABS1 = abs(Y1);
    ABS2 = abs(Y2);

    %sort spectrum from lo to hi to determine the cutoff
    A = sort(abs(Y1));
    B = sort(abs(Y2));
    
    %determine the index that corresponds to cutoff percentage
    cutoffIndex = uint64(length(A)*cutoffPercent);
    
    %determine minimum amplitude to keep in spectrum
    minAmplitude = A(cutoffIndex);

    %eliminate all frequencies below cutoff
    for k = 1:length(A)
        if(ABS1(k) < minAmplitude)
            Y1(k) = gain*Y1(k);
        end
        if(ABS2(k) < minAmplitude)
            Y2(k) = gain*Y2(k);
        end
    end

    %output signal
    z(:,1) = ifft(Y1);
    z(:,2) = ifft(Y2);
    
    %plot, remove comments to enable
    %subplot(2,1,2)
    %stem(abs(fft(z(:,1))))
    %set(gca,'yscal','log')
    %title('output')
end

function z = decon_function(y, Fs, convState)

%     [y,Fs] = audioread(nameIn);
%     samples = [1, len*Fs];
%     clear y Fs;
%     [y,Fs] = audioread(nameIn,samples);

    Fn = Fs/2;                                              % Nyquist Frequency (Hz)
    Wp = 1000/Fn;                                           % Passband Frequency (Normalised)
    Ws = 1010/Fn;                                           % Stopband Frequency (Normalised)
    Rp =    60;                                               % Passband Ripple (dB)
    Rs =    160;                                               % Stopband Ripple (dB)
    [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
    [z,p,k] = cheby2(n,Rs,Ws,'high');                        % Filter Design
    [soslp,glp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability
    % figure(1);
%     freqz(soslp, 2^16, Fs);
%     title('High Pass Filter');
    if(convState == 1)
        disp('Starting convolution handling');
%         figure(2);
%         hold on;
        fd = fft(y);
        li = log(fd);
        outpFreq = filtfilt(soslp,glp,li);
        outpFreqInLog = exp(outpFreq);
        outpFreqFin = ifft(outpFreqInLog);
%         plot(1:length(outpFreq), -outpFreq, 'y');
%         plot(1:length(outpFreqInLog), -outpFreqInLog, 'b');
%         hold off;
        
%         xlim([0 length(outpFreq)]);
%         ylim([0 8]);
%         title('Audio Signal');
%         xlabel('Samples');
%         ylabel('Amplitude');
%         legend({'Original Signal','Noise','Filtered Signal'},'Location', 'northeast');
        x = real(double(outpFreqFin));

    elseif(convState == 2)
        disp('No convolution');
        fd = y;
        li = log(fd);
        outpFreq = filtfilt(soslp,glp,li);
        outpFreqFin = exp(outpFreq);
        x = real(double(outpFreqFin));
    end
    z = y - x;
%     plot(1:length(z), z);
% 	sound(y, Fs);
% 	pause(32);
% 	sound(z, Fs);
    
%     High Pass FIR Filter (Didn't end up using)
%     d = designfilt('highpassfir', ...       % Response type
%        'StopbandFrequency',200, ...     % Frequency constraints
%        'PassbandFrequency',350, ...
%        'StopbandAttenuation',60, ...    % Magnitude constraints
%        'PassbandRipple',4, ...
%        'DesignMethod','kaiserwin', ...  % Design method
%        'ScalePassband',false, ...       % Design method options
%        'SampleRate',2000);
%    fvtool(d);
%     at = filter(hpf,[li; zeros(D,1)]);
%     input = at(D+1:end);
%     z = filter(hpf,input);
    
end