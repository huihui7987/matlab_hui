clear all
R=50e-6;
lamda=1438.5e-9:1e-12:1440.5e-9;
v=(3e8./lamda)-(3e8./1439.5e-9);
%neff=3.17995707;
neff=3.17997;
r=0.83;

FWHM=50e-12;            %��˹�ź�FWHM��ȣ�Ϊ50ps
time_window=100*FWHM;   %��˹�źŵĲ������ڿ�ȣ���ֵ�����˸���Ҷ�任���Ƶ�ʷֱ���
Ns=2001;                %������
dt=time_window/(Ns-1);  %����ʱ����
%dt=1e-12;
t=0:dt:time_window;     %����ʱ��

gauss_time=exp(-0.5*(2*sqrt(2*log(2))*(t-2.5e-9)/FWHM).^2); %��˹���壬����λ��2.5ns����
plot(t*1e+9,gauss_time,'linewidth',2.5);
xlabel('Time/ns');
ylabel('Amplitude/V');
title('Gauss pulse');
%===========���¼���˫���ס�˫�߹����ס�˫�߹������ܶ�=================
%gauss_spec=fftshift(fft(ifftshift(gauss_time)));    %����Ҷ�任�����ҽ�����λ������
gauss_spec=fftshift(fft(ifftshift(gauss_time)));

%%%%0602,test%%%%%%%%%%%%%
%gauss_spec_v = fft(fftshift(gauss_time));֤������0629
gauss_spec_v=fftshift(fft(ifftshift(gauss_time)));
gauss_spec_v=gauss_spec_v/Ns;
%%%%%%%%%%%%%%%%%%%%%%%%%

gauss_spec=gauss_spec/Ns;   %��ʵ�ʵķ���ֵ��

df=1/time_window;               %Ƶ�ʷֱ���
k=floor(-(Ns-1)/2:(Ns-1)/2);    
% k=0:Ns-1;
double_f=k*df;   %˫��Ƶ�׶�Ӧ��Ƶ��,Ƶ�������



L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
taoa=0.83;
Ta1= exp(1i*(pi+phi)).*(taoa-r.*exp(-1i*phi))./(1-r.*taoa.*exp(1i*phi));
%Ta1= (taoa-r)./(1-r.*taoa);
T1= (abs(Ta1)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI1 = angle(Ta1);%~~~~-
if r<=taoa
    PHI1 = PHI1+(PHI1<0)*2*pi ;
end
%PHI=pi+phi+atan((r.*sin(phi))./(tao-r.*cos(phi)))+atan((r*tao.*sin(phi))./(1-tao*r.*cos(phi)))
taob=0.801;
Ta2= exp(1i*(pi+phi)).*(taob-r.*exp(-1i*phi))./(1-r.*taob.*exp(1i*phi));
T2= (abs(Ta2)).^2;%~~~~-
PHI2 = angle(Ta2);%~~

taoc=0.846;
Ta3= exp(1i*(pi+phi)).*(taoc-r.*exp(-1i*phi))./(1-r.*taoc.*exp(1i*phi));
T3= (abs(Ta3)).^2;%~~~~-
PHI3 = angle(Ta3);%~~
if r<=taoa
    PHI3 = PHI3+(PHI3<0)*2*pi ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���亯��,ʵ��n1 = 1.5��
no=1.5;%�˴������޸�n1��ֵ���õ���ͬ����������΢��������
nu=0.4;
nc=1;
H_idea=(1i*2*pi*double_f*1e-9).^(no);

%%%�ֶΣ���ʱ��֪�����������Ŀǰȷ�����ԣ�0512
%H1=(abs(2*pi*double_f*1e-9)).*(n1)*exp(1i*n1*pi/2).*(double_f*1e-9>0)+(abs(2*pi*double_f*1e-9)).*(n1)*exp(-1i*n1*pi/2).*(double_f*1e-9<0);
idea_Intensity_Transmission = mapminmax((abs(1i*2*pi*double_f*1e-9)).^(2*no),0,1);
% idea_power_spec = (abs(H1)).^2;
 idea_phase_spec = angle(H_idea);
 
%%%%%��ͼ%%%%%%
%%%%%%%%%%%%%%%
 %���봫�亯��
%  figure;
%  subplot(1,3,1)
%  plot(double_f*1e-9,abs(H_idea),'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_tran_spec');title('n1�����봫�亯��');hold on;
%  %���봫�亯����λ
%  subplot(1,3,2)
%  plot(double_f*1e-9,idea_phase_spec,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('idea_phase_spec');title('������λ�仯');hold on;
%  %���봫�亯��͸����
% subplot(1,3,3)
% plot(double_f*1e-9,idea_Intensity_Transmission,'r','linewidth',2.5); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');title('����͸����');hold on;


%%%%%��˹��������΢��֮��Ĺ�����%%%%%
% figure;
%����΢�ֽ��
idea_gauss_diff=gauss_spec.*H_idea;
%����΢�ֹ�����
idea_gauss_diff_power_spec=(abs(idea_gauss_diff)).^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ffff =Ta1 .* gauss_spec;


ff_ide  = (1i*2*pi*double_f*1e-9).^(1) .*gauss_spec;
%t0=-3*1e-2;
t0=-0.01;
Hro= -exp(-1i*double_f*1e-9*t0).*(1i*2*pi*double_f*1e-9).^(1);

ff_mo = Hro.* gauss_spec;
% subplot(1,2,1);
% plot(double_f,T1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,T2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% plot(v,T3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Intensity Transmission');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')
% subplot(1,2,2);
% plot(v,PHI1,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% plot(v,PHI2,'b','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% plot(v,PHI3,'g','linewidth',2); xlabel('Frequency(Hz��');ylabel('Phase Response');hold on;
% legend('Critical-coupled','Under-coupled','Over-coupled')
%subplot(1,3,3);plot(v,y3,'linewidth',2);title('n=1.5');
% figure;
% subplot(1,2,1);
% plot(double_f,abs(ffff).^2/3.5/1e-6/0.9397,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('res');hold on;
% subplot(1,2,2);
% plot(double_f,abs(ff_ide).^2,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('res');hold on;

% figure;
% subplot(2,1,1);
% plot(v,Ta2,'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('Ta1');hold on;
% subplot(2,1,2);
% plot(v,abs(gauss_spec),'r','linewidth',2); xlabel('Frequency(Hz��');ylabel('gauss_spec');hold on;


%%%%%%%%%%%%%Ƶ��2ʱ��%%%%%%
x_out_bas = ifftshift(ifft(fftshift(ffff)));%~~~~1��΢��΢�ֽ��ʱ�����в�ͬ
% x_out_bas_vv = ifft(ffff_v);

x_f_ide = ifftshift(ifft(fftshift(ff_ide)));%1������
% x_f_ide_v = ifftshift(ff_ide_v);
x_out_mo = ifftshift(ifft(fftshift(ff_mo)));

figure;

% plot(t*1e+9,abs(x_out_bas)/0.000107,'r','linewidth',2); xlabel('Time/ns');ylabel('Amplitude(a.u.)');hold on;
% plot(t*1e+9,abs(x_f_ide)/0.001852,'g','linewidth',2); xlabel('Time/ns');ylabel('Amplitude(a.u.)');hold on;
%%0.4
% plot(t*1e+9,abs(x_out_bas).^2/(1.145e-8),'r','linewidth',2); xlabel('Time/ns');ylabel('power(a.u.)');hold on;
% plot(t*1e+9,abs(x_f_ide).^2/(3.429e-6),'b','linewidth',2); xlabel('Time/ns');ylabel('power(a.u.)');hold on;
%%1.5
%plot(t*1e+9,abs(x_out_bas).^2/(1.093e-8)/1.048,'k','linewidth',2); xlabel('Time/ns');ylabel('power(a.u.)');hold on;
plot(t*1e+9,abs(x_out_mo).^2/0.0173/(1.972*1e-4)/58.58/1.02,'k','linewidth',2); xlabel('Time/ns');ylabel('power(a.u.)');hold on;
plot(t*1e+9,abs(x_f_ide).^2/0.016/(2.143*1e-4)/59.08,'r','linewidth',2); xlabel('Time/ns');ylabel('power(a.u.)');hold on;

legend('Ring based','ideal')
% figure;
% % 
% plot(t*1e+9,abs(x_f_ide_v),'r','linewidth',2); xlabel('Time/ns');ylabel('Amplitude/V');title('1������΢��ʱ��_v');hold on;
figure;
mmmm=diff(gauss_time);
plot(abs(mmmm),'r','linewidth',2); xlabel('Time/ns');ylabel('Amplitude/V');title('idel');hold on;
plot(abs(mmmm),'b','linewidth',2); xlabel('Time/ns');ylabel('Amplitude/V');title('Ring based');hold on;

