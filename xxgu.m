%�������ڻ��Ƹ�˹������׵������塢���崮��ʱ��Ƶ��ͼ
clc;
clear all;
PlusTime=1;   %����ʱ��ns
PlusType=0;   %��������0:gaus��1:��˹1�׵���2:��˹2�׵�;5��˹5�׵�
PlusSamples=32; %�������������
 
pl=gausplus(PlusSamples,PlusType);
Pn=[1,0,0,0]';  %���Ƶ�����
% Pn=2*randint(16,1)-1;  %DSC_UWB��PN����
% Pn=[0,0,0,0,1,0,0,0,...           %TH_UWB��PN����
%     0,1,0,0,0,0,0,0,...
%     0,0,0,0,0,0,1,0,...
%     0,0,0,0,0,0,0,1,...
%     0,0,1,0,0,0,0,0]';
 
C=pl*Pn';   Tx=C(:);
 
     y=Tx;  TimeRatio=PlusTime/PlusSamples;  %5/7;
     y=y/max(y);  subplot(121);   grid;  plot((1:length(y))*TimeRatio,y,'r-'); xlabel('t(ns)'); legend('ʱ�ǲ���',0); grid on;
     %��y����Ӧ��fft�任,��ͼƵ��ʱ������ȡ����ֵ
     MaxSample=2^nextpow2(length(Tx));  MaxTime=MaxSample*TimeRatio;  fmin=1/MaxTime;
     yfft=fft(y,MaxSample); yfabs=abs(yfft);
     %�任yfabsƵ�ʿ�ʼ��0Ƶ,���±��1��ʼ,����fmin
     fend=6;  %��Ҫ���������Ƶ  
     idx=1:floor(fend/fmin); freq=idx*fmin; subplot(122); 
     semilogy(freq,yfabs(idx)./max(yfabs),'b-')  
%      plot(freq,20*log(yfabs(idx)./max(yfabs)),'b-')  
% %      plot(freq,yfabs(idx),'y-')  
%      plot(freq,yfabs(idx)/max(yfabs),'r-'); xlabel('f(GHz)');legend('Ƶ����',0);