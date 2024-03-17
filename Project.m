%%
%Part1
figure;
%1
[Audio, fs] = audioread("greensleeves.wav");
AudioX = Audio(:,1); 

plot(AudioX(1:325100,1),'DisplayName','Audio');
%soundsc(AudioX(:,1),fs);
%%
%3 
%Yes the signal can be down sampled before the pitch recognition, since 
%the used sampling frequency (which is 44kHz) is about 30 times higher 
%than the signal's maximum relevant frequency 

%4
x_fft = fft(AudioX(:,1),length(AudioX));
plot(fftshift(abs(x_fft)))
title('Fourier Spectrum')
%ylabel('Magnitude')

%9 notes at 325100 th sample 
[auto_corr_y,lags] = xcorr(AudioX(63740:325100,1));
%subplot(2,1,1),plot(AudioX(1:325100,1)) 
%subplot(2,1,2),
% plot(auto_corr_y(325100:length(auto_corr_y)));
% title('Autocorrelation')

%5
N=144;
%spectrogram(AudioX(63740:325095,1), hann(N), 3*N/4, 4*N, fs, 'yaxis');
%title('Sepctogram of the first 9 tones')


%% Part 2

%1
%Tone 1: 38000 to 64000 (in relation to sample number)
%Tone 2: 64000 to 112638
%tone 3: 112638 to 138904
%tone 4: 138904 to 178222
%tone 5: 178222 to 190148
%tone 6: 190148 to 212707
%tone 7: 212707 to 260890
%tone 8: 260890 to 288443
%tone 9: 288443 to 323090

%2
signal = AudioX(38000:64000,1);
x_fft_1st_sample = fft(signal,length(signal));

L = length(signal);
K = ((fs/L*(0:L-1)))';
K(1:L,2) = (abs(x_fft_1st_sample))';
maximum = max(K(:,2));
[x,y] = find(K==maximum);
freq = K(x(1),1);
disp(freq)


%3

%% Part 3

%1
[ismin,ismax,lags] = segmentTone(AudioX,15000);
Algorithm1 = ToneID1stAlgorithm(lags,ismin,ismax,AudioX,fs);
%Algorithm2 = ToneID2ndAlgorithm(lags,ismin,ismax,AudioX,fs);

%2 
%% Question 3
%Audio signal generation
Ts=1/fs;
time = 0.5;
t=[0:Ts:time];
notes = [Algorithm1(1:9)];
sound = zeros(1,int32(time/Ts+1));
for i=1:length(Algorithm1(1:9))
     sound(i,:) = cos(2.*pi.*notes(i).*t)+cos(4.*pi.*notes(i).*t)+cos(8.*pi.*notes(i).*t)+cos(16.*pi.*notes(i).*t);    
end
sig = reshape(sig',9*length(t),1);
reverb = reverberator("DecayFactor",0.9,"SampleRate",fs);
sigout = reverb(sig);
soundsc(sigout(:,1),fs)

%%
%Synthesized signal
Ts=1/fs;
time = 0.8;
t=[0:Ts:time];
notes = [87.31 98 110 98]; 
sound = zeros(1,int32(time/Ts+1));
for i=1:length(notes)
    
        sound(i,:) = cos(2.*pi.*notes(i).*t);
   
end
reverb = reverberator("DecayFactor",0.0,"SampleRate",fs);
sig = reshape(sound',length(notes)*length(t),1);
sigout = reverb(sig);
sigout = sigout./ max(abs(sigout(:)));

%soundsc(sigout(:,1),fs)

% sigout = sigout./ max(abs(sigout(:)));
% audiowrite(filename,sigout(:,1),fs);

[ismin,ismax,lags] = segmentTone(sigout(:,1),15000);
Algorithm1 = ToneID1stAlgorithm(lags,ismin,ismax,sigout(:,1),fs);

%% Functions
%1
function [ismin,ismax,lags] = segmentTone(Signal,Filterwidth)
    figure;
    %Calculates the power of the signal
    PowAudioX = Signal.^2;
    
    %hann window so we calculate the average of the points, giving more importance to the main point we're analysing (remove noise from power signal)
    [mean_pow,lags] = xcorr(PowAudioX,hann(Filterwidth));
    %As the hannis window isn't centered in 0, we need a offset on lags to center it
    lags(:) = lags(:)+(Filterwidth/2);
    
    %finds the index where lags = 0
    zero_lag = find(~lags);
    
    %Find the local minimum
    ismin = islocalmin(round(mean_pow));
    ismin(1:zero_lag)=false;
    
    %Find the local maximum
    ismax = islocalmax(round(mean_pow));
    ismax(1:zero_lag)=false;
    
    
    plot(max(mean_pow).*Signal);
    hold on
    plot(lags(zero_lag:end),mean_pow(zero_lag:end),lags(ismin),mean_pow(ismin),'r*',lags(ismax),mean_pow(ismax),'b*');

end

%% Algorithm for the segmentation of tones and pitches
function f_array = ToneID1stAlgorithm(lags,ismin,ismax,AudioX,fs)
    figure;
    k = find(ismin==1);
    mins = lags(k);
    m = find(ismax==1);
    maxs = lags(m);
    f_array = [0];
    
    for i=1:length(maxs)
        if i>length(mins)
            sig = AudioX(maxs(i):end,1);
        else
            sig = AudioX(maxs(i):mins(i),1);
        end
        fft_note = fft(sig,length(sig));
        L = length(sig);
        
        G = ((fs/L*(0:L-1)))';
        G(1:L,2) = (abs(fft_note))';
       
        maximum = max(G(:,2));
        [x,y] = find(G(:,2)==maximum);
        f_array(i) = G(x(1),1);
    end
    plot((fs/L*(0:L-1)),abs(fft_note));
    fprintf('The frequencies of the recognized tones are, using the first method: \n');
    disp(f_array)
end


%% 2nd Algorithm 
function f_array = ToneID2ndAlgorithm(lags,ismin,ismax,AudioX,fs)

    k = find(ismin==1);
    mins = lags(k);
    f_array = [0];
    
    for i=1:length(mins)-1
        sig = AudioX(mins(i):mins(i+1),1);
    
        fft_note = fft(sig,length(sig));
        L = length(sig);
        
        G = ((fs/L*(0:L-1)))';
        G(1:L,2) = (abs(fft_note))';
       
        maximum = max(G(:,2));
        [x,y] = find(G(:,2)==maximum);
        f_array(i) = G(x(1),1);
    end
    
    fprintf('The frequencies of the recognized tones are, using the second method: \n')
    disp(f_array)
end

%% 2
function Error_array = ErrorCalculator(GroundTruth,Freq_array)
    Error_array = abs(Freq_array-GroundTruth);
    fprintf('Mean error: \n')
    disp(mean(Error_array))
end
%%GroundTruth = [329.63,392,440,493.23,523.25,493.88,440,369.99];


