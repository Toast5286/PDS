%%
%Part1

[Audio, fs] = audioread("greensleeves.wav");
AudioX = Audio(:,1); 

plot(AudioX,'DisplayName','Audio');
%soundsc(AudioX(:,1),fs);

%3 
%Yes, the signal can be down sampled before the pitch recognition, since
%the it's maximum relevant frequency is about 30 times higher than the used
%sampling frequency (which is 44kHz)


%4
x_fft = fft(AudioX(:,1),length(AudioX));
%plot(fftshift(abs(x_fft)))

%9 notes at 325095 th sample 
[c,lags] = xcorr(AudioX(176880:188574,1),AudioX(1:325095,1));
stem(lags,c); 

%5
N=144;
%spectrogram(AudioX(1:325095,1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
%title('Sepctogram of the first 9 tones')

%% Part 2


%% Part 3


