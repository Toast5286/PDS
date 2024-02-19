%% R1.a
sRate = 8000;
sInterval = 1/sRate;
t1 = 0:sInterval:2;
[xSampled,instFreq] = xc(t1);
plot(t1,xSampled,'ro-')

maxInstFreq = 2*pi*(1000*(2).^2);
fprintf('The maximum instant frequency is: %d', maxInstFreq)

soundsc(x, 8000); 

%Comment 

%% R2.a

x = xSampled;
N=72;
spectrogram(x, hann(N), 3*N/4, 4*N, 8000, 'yaxis');

%Comment
%The spectogram indicates that with time the frequency of the signal
%becomes higher, accordingly the perceived pitch of the signal/sound also
%gets higher 

%% R3.a

%Since the sampling frequency of x(n) is 8 kHz, (8000 samples per second),
%then the sampling frequency of y(n) is 4 kHz (4000 samples per second).  

sRate2 = 4000;
sInterval2 = 1/sRate;
t2 = 0:sInterval2:2;
[xSampled2,instFreq2] = xc(t1);
plot(t2,xSampled2,'ro-')

maxInstFreq = 2*pi*(1000*(2).^2);
fprintf('The maximum instant frequency is: %d', maxInstFreq)

N=72;
spectrogram(xSampled2, hann(N), 3*N/4, 4*N, 4000, 'yaxis');
soundsc(xSampled2, 4000);  

% R3.b

sRate_z1= 20000;
sInterval_z1 = 1/sRate_z1;
tz1 = 0:sInterval_z1:2;
[xSampled_z1,instFreq2] = z1(tz1);

sRate_z2= 8000;
sInterval_z2 = 1/sRate_z2;
tz2 = 0:sInterval_z2:2;
[xSampled_z2,instFreq2] = z2(tz2);

sRate_z3= 4000;
sInterval_z3 = 1/sRate_z3;
tz2 = 0:sInterval_z3:2;
[xSampled_z3,instFreq2] = z3(tz3);



%Comment 

% R3.c




%% R4.a

%% R5.a

%% R6.a



%% Function Definitions
function x = z1(t)
f1 = 5000;
x = cos(2*pi*f1*t);
end
function x = z2(t)
f2 = 2500;
x = cos(2*pi*f2*t);
end
function x = z3(t)
f3 = 1000;
x = cos(2*pi*f3*t);
end







