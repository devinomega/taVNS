intensity = 1;
pulseW = .0005;     %Pulse width in seconds
sRate = 50000;      %Sampling frequency

% Burst
repFreq = 5;        
burstFreq = 100;   %Pulse frequency in Hz %Burst Freq
numBurst = 5; 
t = 0 : 1/sRate : 1/repFreq; % sampling frequency 
d = 0: 1/burstFreq : (numBurst-1)/burstFreq; %repetition frequency   
y =intensity*pulstran(t,d,'rectpuls',pulseW);
y(end) = [];
wave_burst = repmat(y',repFreq,1);

%25 Hz
repFreq = 25;       %Pulse frequency in Hz
t = 0 : 1/sRate : 1; % sampling frequency   
d = pulseW/2: 1/repFreq : 1; %repetition frequency   
y_25 =intensity*pulstran(t,d,'rectpuls',pulseW);
y_25(end)= [];

   
%100 Hz
repFreq = 100;      %Pulse frequency in Hz
t = 0 : 1/sRate : 1; % sampling frequency   
d = pulseW/2: 1/repFreq : 1; %repetition frequency   
y_100 =intensity*pulstran(t,d,'rectpuls',pulseW);
y_100(end)= [];

        
t(end) = [];
t = t'*1000;

figure 
hold on

subplot(3,1,1)
plot(t,wave_burst,'LineWidth',1.5)
ylim([0 1.5])
xlabel('Time (ms)')
ylabel('Intensity (mA)')
title('Burst Waveform')
box off

subplot(3,1,2)
plot(t,y_25,'LineWidth',1.5)
ylim([0 1.5])
xlabel('Time (ms)')
ylabel('Intensity (mA)')
title('25 Hz Waveform')
box off

subplot(3,1,3)
plot(t,y_100,'LineWidth',1.5)
ylim([0 1.5])
xlabel('Time (ms)')
ylabel('Intensity (mA)')
title('100 Hz Waveform')
box off
