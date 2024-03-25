%%
%Project done by:
    % Gil Beirão nº99943
    % Rodolfo Amorim nº100260

%Index
%Part 1 from line 11 to line 74 
%Part 2 from line 75 to line 117 
%Part 3 from line 118 to line 194
%Functions section from line 195 to 301
%% Part1
figure;
%1
[Audio, fs] = audioread("greensleeves.wav");
AudioX = Audio(:,1); 
%soundsc(AudioX(:,1),fs);

plot(AudioX);
title('\textbf{Signal from greensleeves.wav}', 'Interpreter','latex')
xlabel('\textbf{Sample}','Interpreter','latex');
ylabel('\textbf{Amplitude}', 'Interpreter','latex');

%2
plot(AudioX(1:325100,1));
title('\textbf{Shape of first 9 tones}', 'Interpreter','latex')
xlabel('\textbf{Sample}','Interpreter','latex');
ylabel('\textbf{Amplitude}', 'Interpreter','latex') ;

%3 
%Yes, the signal can be down sampled before the pitch recognition, since 
%the used sampling frequency (which is 44KHz) is about 10 times higher 
%than the signal's maximum relevant frequency 

%4
%Tone 1: 38000 to 64000 (in relation to sample number)
%Tone 2: 64000 to 112638
%tone 3: 112638 to 138904
%tone 4: 138904 to 178222
%tone 5: 178222 to 190148
%tone 6: 190148 to 212707
%tone 7: 212707 to 260890
%tone 8: 260890 to 288443
%tone 9: 288443 to 323090

list1 = [38000 64000 112638 138904 178222 190148 212707 260890 288443]; %beginning of each note
list2 = [64000 112638 138904 178222 190148 212707 260890 288443 325100]; %ending of each note 

%Fourier Spectrum
for i=1:length(list1)
    signal = AudioX(list1(i):list2(i),1);
    L = length(signal);
    sig_fft = fft(signal,length(signal));
    axis = ((fs/L*(0:L-1)));
    subplot(3, 3, i);
    plot(axis,abs(sig_fft));
    title('\textbf{Fourier Spectrum of note }' + string(i), 'Interpreter','latex')
    xlabel('\textbf{Frequency (Hz)}','Interpreter','latex');
    ylabel('\textbf{Magnitude}', 'Interpreter','latex');
end

%Auto-Correlation
for r=1:length(list1)
    [auto_corr,lag] = xcorr(AudioX(1:325100,1),AudioX(list1(r):list2(r),1));
    subplot(3, 3, r);
    plot(auto_corr(325100:length(auto_corr)));
    title('\textbf{Auto-Correlation of note }' + string(r), 'Interpreter','latex')
    xlabel('\textbf{Lags}','Interpreter','latex');
    ylim([0 1000]);
end

%5
N=144;
spectrogram(AudioX, hann(N), 3*N/4, 4*N, fs, 'yaxis');
title('Sepctogram of the first 9 tones')
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
% clc
% signal = AudioX(38000:64000,1);
% x_fft_1st_sample = fft(signal,length(signal));
% L = length(signal);
% K = ((fs/L*(0:L-1)))';
% K(1:L,2) = (abs(x_fft_1st_sample))';
% maximum = max(K(:,2));
% [x,y] = find(K==maximum);
% freq = K(x(1),1);
% fprintf('The pitch recognition value is = %f \n',freq);

%3
clc
Ts=1/fs;
time = 0.7;
t=[0:Ts:time];
notes = [440];
signal = zeros(1,int32(time/Ts+1));
for i=1:length(notes)
        signal(i,:) = cos(2.*pi.*notes(i).*t)+ cos(16.*pi.*notes(i).*t)+ cos(32.*pi.*notes(i).*t);
end
sig_fft = fft(signal,length(signal));
L = length(signal);
K = ((fs/L*(0:L-1)))';
K(1:L,2) = (abs(sig_fft))';
maximum = max(K(:,2));
[x,y] = find(K==maximum);
freq = K(x(1),1);
fprintf('The absolute error for the pitch recognition is = %f \n',notes(1)-freq);
%% Part 3

%1
[ismin,ismax,lags] = segmentTone(AudioX(:,1),fs,true);
Notes_syn = ToneID1stAlgorithm(lags,ismin,ismax,AudioX,fs);

%2
% Segmentation algorithm Evalutaion 
TimeIntervale = 0.2:0.2:1;
Decay_Factor = [100 80 60 40 20 10 5 2];
AbsoluteError = zeros(length(TimeIntervale),length(Decay_Factor));
freq = 200;
EndTime = 100;

for TimeIntervale_index = 1:length(TimeIntervale)
    Starts = 1:TimeIntervale(TimeIntervale_index):EndTime-TimeIntervale(TimeIntervale_index);
    Ends = 1+TimeIntervale(TimeIntervale_index):TimeIntervale(TimeIntervale_index):EndTime;
    
    for decay_index = 1:length(Decay_Factor)
        sound = zeros(1,ceil(Ends(end)*fs));
        
        for time_index = 1:length(Starts)
            t = 0:(1/fs):(Ends(time_index)-Starts(time_index));
            SampleStart = round(Starts(time_index)*fs);
            
            sound(1,SampleStart+1:(SampleStart+length(t))) = cos(2*pi*freq*t).*exp(-Decay_Factor(decay_index)*t);
        end
    
        [ismin,ismax,lags] = segmentTone(sound,fs,false);
        ismin(end)=1;
        Predicted_Starts = lags(ismax)/fs;
        Predicted_Ends = lags(ismin)/fs;
        if length(Starts) == length(Predicted_Starts) & length(Ends) == length(Predicted_Ends)
            AbsoluteError(TimeIntervale_index,decay_index) = mean([ErrorCalculator(Starts,Predicted_Starts),ErrorCalculator(Ends,Predicted_Ends)])/TimeIntervale(TimeIntervale_index);
            
        end
    end
end

plot(Decay_Factor,AbsoluteError);
legend('Duration: 0.2s','Duration: 0.4s','Duration: 0.6s','Duration: 0.8s','Duration: 1.0s' );
xlabel("Decay Factor (Hz)");
ylabel("Relative Mean Absolute Error (%)");

%General Evaluation
%Evaluation for the first 9 notes of the greensleeves audio file
% E G A B C B A F# D

plot(AudioX(1:323090,1));
GroundTruth_Freq = [329.63 392 440 483.88 523.25 493.88 440 369.99 293.66];

GroundTruth_Start = [38000 64000 112638 138904 178222 190148 212707 260890 288443];
GroundTruth_End = [64000 112638 138904 178222 190148 212707 260890 288443 323090];

[ismin,ismax,lags] = segmentTone(AudioX(1:323090,1),fs,true);
Algorithm1 = ToneID1stAlgorithm(lags,ismin,ismax,AudioX,fs);

Pred_Start = lags(ismax);
Pred_End = lags(ismin);

FreqErrorArray = ErrorCalculator(GroundTruth_Freq,Algorithm1);
StartErrorArray = ErrorCalculator(GroundTruth_Start,Pred_Start)/fs;
EndErrorArray = ErrorCalculator(GroundTruth_End,Pred_End)/fs;

%% 
%Synthesised Greensleeves signal with the fundamental frequencies 
%given by the previous Algorithm 
Ts=1/fs;
time = 0.59;
t=[0:Ts:time];
notes = Notes_syn;
sound = zeros(1,int32(time/Ts+1));
for i=1:length(notes)
        sound(i,:) = ((cos(2.*pi.*notes(i).*t)+cos(4.*pi.*notes(i).*t))).*exp(-5.*t);
end
sig = reshape(sound',length(notes)*length(t),1);
soundsc(sig,fs);
%% Functions

%Segmentation Algorithm
function [ismin,ismax,lags] = segmentTone(Signal,fs,plot_plz)
    Filterwidth = round(fs*15000/44000);

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
    ismin(zero_lag+length(Signal):end)=false;
    ismin(zero_lag+length(Signal)) = true;
    
    %Find the local maximum
    ismax = islocalmax(round(mean_pow));
    ismax(1:zero_lag)=false;
    ismax(zero_lag+length(Signal):end)=false;
    
    if(plot_plz)
        figure;
        plot(max(mean_pow).*Signal);
        hold on
        plot(lags(zero_lag:end),mean_pow(zero_lag:end),lags(ismin),mean_pow(ismin),'r*',lags(ismax),mean_pow(ismax),'b*');
    end
end

%Algorithm for the recognition of pitches using DFT
function f_array = ToneID1stAlgorithm(lags,ismin,ismax,AudioX,fs)
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

    figure;
    plot((fs/L*(0:L-1)),abs(fft_note));
    fprintf('The frequencies of the recognized tones are, using the first method: \n');
    disp(f_array)
end

%Error Calculator for the absolute error metric 
function Error_array = ErrorCalculator(GroundTruth,Freq_array)
    Error_array = abs(Freq_array-GroundTruth);
    fprintf('Mean error: \n')
    disp(mean(Error_array))
end

%Algorithm for the recognition of pitches using Auto-Correlation
function f_array = ToneID3rdAlgorithm(lags,ismin,ismax,AudioX,fs)
    figure;
    k = find(ismin==1);
    mins = lags(k);
    m = find(ismax==1);
    maxs = lags(m);
    f_array = [0];
    sig_Period_sample = [0];
    LMax_index = [0,0];
    
    %Samples/S * 1/s (smallest/biggest frequency a human can hear)
    minSampleLag = fs/19.44;
    
    for i=1:length(maxs)
        if i>length(mins)
            sig = AudioX(maxs(i):end,1);
        else
            sig = AudioX(maxs(i):mins(i),1);
        end
        [auto_corr,auto_lag] = xcorr(sig,sig,ceil(minSampleLag));
        LocalMax = auto_corr(islocalmax(auto_corr));
        LMax_index(1) = min(find(auto_corr == LocalMax(1)));
        LMax_index(2) = min(find(auto_corr == LocalMax(2)));
        sig_Period_sample(i) = LMax_index(2)-LMax_index(1);

        f_array(i) = (fs/sig_Period_sample(i));
    end
    
    plot(auto_lag,auto_corr);
    fprintf('The frequencies of the recognized tones are, using the third method: \n');
    disp(f_array)
end