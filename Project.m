%%
%Part1

[Audio, fs] = audioread("greensleeves.wav");
AudioX = Audio(:,1); 

x_fft = fft(AudioX(:,1),length(AudioX));
plot(fftshift(abs(x_fft)))

plot(AudioX,'DisplayName','Audio');
soundsc(AudioX(:,1),fs)

%%