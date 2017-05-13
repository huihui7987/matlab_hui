clear all
FWHM = 10;
row = FWHM/(2*sqrt(2*2*log(2)));
begin = 1430;
stop = 1450;
center = (stop+begin)/2;
for i = begin:stop
wavelength(i-begin+1)=i;
fun(i-begin+1)=exp(-(i-center)^2/(2*row*row));
end
plot(wavelength,fun)
%%

dt=.0001;
fs=1/dt; %sampling frequency
fn=fs/2;
n=1000;
t=dt*(-n/2:n/2); %time base

sigma=0.001;
variance=sigma^2;

f0 = 1000;
xt=cos(2*pi*t*f0) .* (exp(-t.^2/(2*variance)))/sqrt(2*pi*variance);
subplot(2,1,1); plot(t,xt,'b'); 
title(['Gaussian Pulse \sigma=', num2str(sigma),'s']);
xlabel('Time(s)'); ylabel('Amplitude');
axis([-0.02 0.02 -400 400]);

xf = fftshift(fft(xt));
f = fs*(-n/2:n/2)/(n/2); %Frequency Vector
subplot(2,1,2); plot(f,abs(xf),'r'); title('Magnitude of FFT');      
xlabel('Frequency (Hz)'); ylabel('Magnitude |X(f)|');