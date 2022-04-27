% Lowpass filter for 0-20 Hz
lp20 = designfilt('lowpassfir', 'PassbandFrequency', 20, 'StopbandFrequency', 30,...
    'PassbandRipple', 0.5, 'StopbandAttenuation', 60, 'SampleRate', 70);

% Bandpass filter for 20-60 Hz
bp20_60 = designfilt('bandpassfir','FilterOrder',60,'CutoffFrequency1',20,...
   'CutoffFrequency2',60,'StopbandAttenuation1',60,'StopbandAttenuation2',60,'SampleRate',200);

% Bandpass filter for 60-200 Hz
bp60_200 = designfilt('bandpassfir','FilterOrder',69,'CutoffFrequency1',60,...
    'CutoffFrequency2',200,'StopbandAttenuation1',60,'StopbandAttenuation2',60,'SampleRate',500);

% Bandpass filter for 200-600 Hz
bp200_600 = designfilt('bandpassfir','FilterOrder',120,'CutoffFrequency1',200,...
    'CutoffFrequency2',600,'StopbandAttenuation1',60,'StopbandAttenuation2',60,'SampleRate',1300);

% Bandpass filter for 600-3000 Hz
bp600_3k = designfilt('bandpassfir','FilterOrder',167,'CutoffFrequency1',600,...
    'CutoffFrequency2',3000,'StopbandAttenuation1',60,'StopbandAttenuation2',60,'SampleRate',6100);

% Bandpass filter for 3000-8000 Hz
bp3k_8k = designfilt('bandpassfir','FilterOrder',330,'CutoffFrequency1',3000,...
    'CutoffFrequency2',6000,'StopbandAttenuation1',60,'StopbandAttenuation2',60,'SampleRate',12100);

% Highpass filter for >8000 Hz
hp8k = designfilt('highpassfir','CutoffFrequency',8000,'FilterOrder',440,...
     'PassbandRipple',0.5,'StopbandAttenuation',60,'SampleRate', 16100);

% All code below is for testing purposes, not to be used 
% lp = designfilt('lowpassfir', 'PassbandFrequency', 225, 'StopbandFrequency', 275,...
%     'PassbandRipple', 0.5, 'StopbandAttenuation', 200, 'SampleRate', 2000);
% bp = designfilt('bandpassfir','FilterOrder',100,'CutoffFrequency1',250,...
%     'CutoffFrequency2',2000,'SampleRate',5000);
% hp = designfilt('highpassfir','StopbandFrequency',2000,'PassbandFrequency',2050,...
%     'PassbandRipple',0.5,'StopbandAttenuation',65,'SampleRate', 5000);

% fvtool(hp8k)  
