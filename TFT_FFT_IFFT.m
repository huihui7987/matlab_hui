% temperal differiator
close all;
clear;
clc;
c=3e8; 

%%%%%%%%% 
N = 1000+1; 
t_length = 10;

t = linspace(0,t_length,N);
t_interval = t_length/(N-1);

%%%%%%%%% 
wavelength_carrier = 1550e-9;
f_carrier = c/wavelength_carrier;
f_interval = 1/t_length; 

f_bas = -(N-1)/2*f_interval:f_interval:(N-1)/2*f_interval; %Ƶ������

%%%%%%%%% 
%x_in_bas = sin(2*pi*t*0.2); 
x_in_bas = sin(t);% 
%x_in_bas = sin(t).^2.*cos(t).*t;
xF_in_bas = fftshift(fft(x_in_bas)); 

%%%%%%%%% 
w=2*pi*f_bas;%~~~~  
hF = (1i.*w);%~~~~ �е���΢�������ݺ���
%hF = 1;% 
%hF = exp(-1i*w*2.5);% 
xF_out_bas = xF_in_bas.*hF;%~~~~ %΢�����

%%%%%%%%% 
x_out_bas = ifft(ifftshift(xF_out_bas));%~~~~ ΢�ֽ��ʱ��

%%%%%%%%% 
x_out_diff_test = [diff(x_in_bas), 0]./t_interval;


figure(1)
subplot(3,2,1)

plot(t,x_in_bas)
title('����ʱ��');

subplot(3,2,2)
plot(f_bas,abs(xF_in_bas))
title('����Ƶ��');
subplot(3,2,3)
plot(f_bas,abs(xF_out_bas))
title('Ƶ�����');
subplot(3,2,4)
plot(t,x_out_bas)
title('ʱ�����');
subplot(3,2,5)
plot(f_bas,abs(xF_in_bas))
title('����Ƶ��');

subplot(3,2,6)
plot(t,x_out_diff_test)
title('����������');

figure(2)
subplot(3,2,1)
plot(f_bas,abs(hF))
title('΢��������');
subplot(3,2,2)
plot(f_bas,angle(hF)./pi)
title('΢������λ');



