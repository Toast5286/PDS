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
[auto_corr_y,lags] = xcorr(AudioX(30000:325100,1));%,AudioX(63740:325100,1));
%auto_corr_y = auto_corr_y(length(AudioX(1:325100,1)):length(AudioX(1:325100,1))+3000);
subplot(2,1,1),plot(AudioX(1:325100,1)) 
subplot(2,1,2),plot(auto_corr_y(325100:length(auto_corr_y)))
%subplot(3,1,3),stem(lags,auto_corr_y) set the subplot(3,...,...) in order
%to show these 3 plots 
%stem(lags,c); 

%5
N=144;
%spectrogram(AudioX(1:325095,1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
%title('Sepctogram of the first 9 tones')


%% Part 2

%1
%Tone 1: 38000 to 64000 (in relation to sample number)

%% Part 3


