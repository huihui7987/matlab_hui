%% Fourier transform of Gaussian pulse
% The code is written by Kh.Tsendsuren, 2015-07-03

clear; close all; clc
c = 3e+8; % Speed of light [m/sec]
lambda = 800e-9; % Wavelength [nm]
freq = c/lambda; % Actual Frequency of light [THz]
%%
% 
%  Nyquist sampling theorem
%  The sampling frequency should be at least twice the highest frequency
%  contained in the signal, fs>2*f
% 

fsamp = freq*10; % Sampling frequency

T = 1/fsamp; % Unit time [fs]
L = 200; % Length of signal
sigma = 8e-15; % Pulse duration
t = (0:L-1)*T; % Time base
t0 = max(t)/2; % Used to centering the pulse

%% Electric field
Egauss = (exp(-2*log(2)*(t-t0).^2/(sigma)^2)).*cos(-2*pi*freq*(t-t0));

subplot(2,1,1)
plot(t/1e-15,real(Egauss),'b');
title(['Gaussian Pulse \sigma=', num2str(sigma),'s']);
xlabel('Time (fs)');
ylabel('Amplitude');
ylim([-1 1]) 
%xlim([30e-15 70e-15])
grid on

NFFT = 2^nextpow2(L);
%NFFT = 1900; % Frequency sampling number, it will define how 
X = fft(Egauss,NFFT)/L;
%Pxx=X.*conj(X)/(NFFT*NFFT); %computing power with proper scaling
freq = 0.5*fsamp*linspace(0,1,NFFT/2+1); % (full range) Frequency Vector
subplot(2,1,2)
plot(freq/1e+12,2*abs(X(1:NFFT/2+1)))
title('Magnitude of FFT');
xlabel('Frequency (THz)')
ylabel('Magnitude |X(f)|');
%xlim([340 380])
grid on