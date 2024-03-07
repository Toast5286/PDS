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
plot(fftshift(abs(x_fft)))

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

%Calculates the power of the signal
PowAudioX = AudioX.^2;

%hann window so we calculate the average of the points, giving more importance to the main point we're analysing (remove noise from power signal)
[mean_pow,lags] = xcorr(PowAudioX(1:325100),hann(15000));
%As the hannis window isn't centered in 0, we need a offset on lags to center it
lags(:) = lags(:)+(15000/2);

%finds the index where lags = 0
zero_lag = find(~lags);

%Find the local minimum
ismin = islocalmin(mean_pow);
ismin(1:zero_lag)=false;

%Find the local maximum
ismax = islocalmax(mean_pow);
ismax(1:zero_lag)=false;


plot(500.*AudioX(1:325100))
hold on
plot(lags(zero_lag:end),mean_pow(zero_lag:end),lags(ismin),mean_pow(ismin),'r*',lags(ismax),mean_pow(ismax),'b*');



% 
% 
% peak_index=[];
% peak_pow = [];
% t=3;
% for i=1:length(mean_pow)-1
%     peak_pow = [peak_pow,mean_pow(i+1)-mean_pow(i)];
%     % if mean(mean_pow(i:i+20))-mean(mean_pow(i-20:i))> t
%     %     peak_index = [peak_index,i];
%     % end
% end
% plot(mean_pow);
% hold on;
% scatter(peak_pow,'r');

