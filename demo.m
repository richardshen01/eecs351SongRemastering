[a, Fs] = demo_functions('sample_one_original.wav',2,100000,1);
a = Second_Gaussian(a,Fs);
a = Sharpening_Filter(a,Fs,'delta_length_1323000.wav');
audiowrite('demo_out.wav',a,Fs);