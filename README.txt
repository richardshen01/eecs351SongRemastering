DEPENDENCIES:
- filter_design.m: initializes filters needed for the equalizer when called.
- test_functions.m: contains functions that are called in the remastering program to apply our dsp techniques.
- Second_Gaussian.m: a function that applies a second Gaussian window after our test_functions.
- Sharpening_Filter.m: a function that uses the impulse response specified in the UI and applies it to the song.
- delta_length_720000.wav: an impulse response .wav file that can be used by Sharpening_Filter.m.
- delta_length_1323000.wav: an impulse response .wav file that can be used by Sharpening_Filter.m.
- sample_one_original.wav: A sample input for our program, but any file can be used.
- demo.m: A short demo of our program that takes sample_one_original.wav and remasters it. It's output is demo_out.wav.
- eq_app.mlapp: the user interface that has the equalizer and settings and is used to generate the remastered song.

HOW TO USE THE APP:
(FIRST ENSURE THAT ALL THE DEPENDENCIES ARE TOGETHER IN THE SAME LOCAL DIRECTORY)
1. Double-click on eq_app.mlapp (MATLAB may open if it is not already open)
2. In the bottom right corner, type in the input file name with the file extension
3. In the "Length" box, type in the desired length of the song you would like to use (this must be less than or equal the actual length of the song)
4. In the "Gain" box, type the desired gain amplification for the overall song. This is set to 1 by default.
5. In the "Window Size" box, type the desired window size to be used. This is set to 100000 by default and is recommended.
6. Select which impulse response you would like to use. For sample_one_original.wav, we recommend using delta_length_1323000.wav.
7. Click "Generate" and the status next to the button should change to "Running."
8. Once the status changes to "Completed," a file called "song_output.wav" should be generated in your local directory.