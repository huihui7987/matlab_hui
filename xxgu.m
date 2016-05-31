%程序用于绘制高斯及其各阶导数脉冲、脉冲串的时、频谱图
clc;
clear all;
PlusTime=1;   %脉冲时间ns
PlusType=0;   %脉冲类型0:gaus；1:高斯1阶导；2:高斯2阶导;5高斯5阶导
PlusSamples=32; %单脉冲抽样次数
 
pl=gausplus(PlusSamples,PlusType);
Pn=[1,0,0,0]';  %绘制单脉冲
% Pn=2*randint(16,1)-1;  %DSC_UWB的PN序列
% Pn=[0,0,0,0,1,0,0,0,...           %TH_UWB的PN序列
%     0,1,0,0,0,0,0,0,...
%     0,0,0,0,0,0,1,0,...
%     0,0,0,0,0,0,0,1,...
%     0,0,1,0,0,0,0,0]';
 
C=pl*Pn';   Tx=C(:);
 
     y=Tx;  TimeRatio=PlusTime/PlusSamples;  %5/7;
     y=y/max(y);  subplot(121);   grid;  plot((1:length(y))*TimeRatio,y,'r-'); xlabel('t(ns)'); legend('时城波形',0); grid on;
     %求y其相应的fft变换,划图频谱时将幅度取绝对值
     MaxSample=2^nextpow2(length(Tx));  MaxTime=MaxSample*TimeRatio;  fmin=1/MaxTime;
     yfft=fft(y,MaxSample); yfabs=abs(yfft);
     %变换yfabs频率开始是0频,而下标从1开始,递增fmin
     fend=6;  %需要划到的最高频  
     idx=1:floor(fend/fmin); freq=idx*fmin; subplot(122); 
     semilogy(freq,yfabs(idx)./max(yfabs),'b-')  
%      plot(freq,20*log(yfabs(idx)./max(yfabs)),'b-')  
% %      plot(freq,yfabs(idx),'y-')  
%      plot(freq,yfabs(idx)/max(yfabs),'r-'); xlabel('f(GHz)');legend('频域波形',0);