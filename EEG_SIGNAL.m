clc; clear

% Display the signal in time and frequency

load("EEG_X");
fs=72.5;
N=length(EEG);

figure("name","EEG Signal & Spectrum","Position",get(0,'ScreenSize'),"Visible","on");
time=(0:(N-1))/fs;
figure(1);
subplot(2,1,1);
plot(time,EEG);
title('EEG Signal');
xlabel('Time [sec]',"FontWeight","bold");
ylabel('Amplitude',"FontWeight","bold");

subplot(2,1,2);
y=abs(fftshift(fft(EEG)));
f_vec=linspace(-fs/2,fs/2-fs/N,N);
stem(f_vec,y);
title('Spectrum of EEG Signal');
xlabel('Frequency [Hz]',"FontWeight","bold");
ylabel('Amplitude',"FontWeight","bold");


%% 
% Calculation of alpha wave frequency in 3 different ways

% Finding the frequency from the time plane

EEG1=EEG-mean(EEG);
%plot(time,EEG1);
%hold on;
dt=1/fs;
[piks, locs]=findpeaks(EEG1,fs,"minpeakdistance",0.05);
%plot(locs,piks,"r*");
diff_picks=mean(diff(locs));
tederA=1/diff_picks;
disp("The frequency is "+ tederA);

% Finding the frequency after a shift to the frequency plane

y1=y((N/2)+2:end);
f1=f_vec((N/2)+2:end);
%stem(f1,y1);
[value, loc]=max(y1);
tederB=f1(loc);
disp("The frequency is "+ tederB);

% Finding the frequency by Z transform

% a - By the equation
n=0:N-1;
% X(z)=sum(EEG.*z.^(-n)) Z transform
dw=2*pi/N;
w=0:dw:2*pi-dw;

X=zeros(1,N);
for i=1:length(w)
    z=exp(-1i*w(i));
    X(i)=sum(EEG.*(z.^n)); % z transform for the unit circle
end

f_vec=linspace(-fs/2,fs/2-fs/N,N);
EEGspec=fftshift(abs(X));
f_vector=f_vec(N/2+2:end);
EEGspectrum=EEGspec(N/2+2:end);
%plot(f_vector,EEGspectrum)

[value1, loc1]=max(EEGspectrum);
tederC=f_vector(loc1);
disp("The frequency is "+ tederC);

% b- by using function freqz

[val,f]=freqz(EEG-mean(EEG),1,2048,fs); 

%figure;
%plot(f,abs(val))
[value2, loc2]=max(val);
tederD=f(loc2);
disp("The frequency is "+ tederD);

%% 
% Butterworth filter on th EEG signal

fc=(0.1*fs)/(0.5*fs); 
[d,c]=butter(2,fc,"high");
EEGfilt=filtfilt(d,c,EEG);

figure("name","Filtered EEG Signal","Position",get(0,'ScreenSize'),"Visible","on");

subplot(2,1,1);
time=(0:(N-1))/fs;
plot(time,EEG);
xlabel('Time [sec]',"FontWeight","bold");
ylabel("EEG Signal","FontWeight","bold");
title("Original EEG Signal");

subplot(2,1,2);
plot(time,EEGfilt);
xlabel('Time [sec]',"FontWeight","bold");
ylabel("EEG Signal","FontWeight","bold");
title("Filtered EEG Signal");








