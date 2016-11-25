%�������ģ������ֵģ������ĵ�
clear all;
clc
c       = 3e8;%����
h       = 6.626e-34;%���ʿ˳���

% TUNING PARAMETERS
P       = 1;                                                               %2; % Watt  [0.2 0.5 1 3 5 10] ���ֹ���       %P=P;%-0.005*P;
BR      = 10e9;                                                            %10e9;          % Bit rate  [1 10 20 40]
PW      =0.5e-12;                                                          %.5e-12;  %1e-12;          % Pulse Width [1 3 5 10 20]
D       =10;%1400;                                                             %1400; %200;-2000-1100;  % Dispersion in ps/nm.km,[-600 -300 100 300 600]
FName   = 'resonatorOL-4.txt';
%%%%%%%%% WAVEGUIDE PARAMETERS
L       =135e-6*2*pi;                                                      %4.4e-5;     %135e-6*2*pi;    %4.4e-4;      % Length in meters
Aeff    = .50e-12;                                                         % Effective Area in m^2
neff    =1.7;                                                              %3.37;  %GaAs   % D=500ps/nm.km ==> B2 = -637.28 ps^2/km
n2      = 6.5e-18;                                                         % (4-9)e-18 m^2/W
B       =0.75e-11;                                                         % TPA (m/W)
tau     =100e-12;%1.0e-12;                                                 %300e-12;     % Free carrier lifetime
alp_dB  = 200;%100;                                                        %2000       % dB/m loss
att     = alp_dB./4.343;

%%%%%%%%% SIMULATION PARAMETERS
I       = P./Aeff;                                                         %���ֹ�ǿ��
BW      = 20;                                                              %�źŴ���
T       = 1/BR;                                                            %��������
Lp      = 1545.55e-9;                                                       % Pump Wavelength
% Ls      = Lp+15e-9;                                                      % Signal Wavelength
% Lc      = Lp-15e-9;                                                      % Conjugate Wavelength
Ls      = Lp + 1e-9*linspace(0.5,BW,1191);                                 % Signal Wavelength
Lc      = Lp - 1e-9*linspace(0.5,BW,1191); %2*Lp-Ls;%                       % Conjugate Wavelength % As = 1e-8;   % % Stokes power    % Ac = 0;  % Conjugate power
As      = 5.13e-5*ones(1,1191); %1e-8*ones(1,1191);                         % Stokes power
Ac      = zeros(1,1191);
wp      = 2*pi*c/Lp;                                                       % angular frequency
ws      = 2*pi*c./Ls;
ga      = n2*wp/(c*Aeff);                                                  %������ϵ��
B2      = -637.28/500*D*1e-27;                                             % s^2/m Ⱥ�ٶ�ϵ��
En      = h*c/Lp;                                                          %��������
As0     = As;
Pc0     = abs(Ac).^2;
Gain    = 1;
N       = 1e-6*(1/(1-exp(-T/tau)))*B*PW/(2*En)*I.^2;                       %1e-6*(tau/T + 0.5)*B*PW/(2*En)*I.^2; �������ܶ�
alfa    = 1.45e-17*N*(Lp./Lp).^2 ;                                          %N = 1e-6*B*PW/(2*En)*I.^2;
loss    = 1.45e-17*N;                                                      %���ֹ��������������
dB      = B2*(wp-ws).^2;                    % Delta K term ��λʧ��
Ns      = 10;
%ȡ������
dz      = L/Ns;%��λȡ������
cAc     = conj(Ac);%Ac�Ĺ���
P1=P;%���ֹ���

Ngain=1e-20*ones(1,1191);
Nloss=1e-20*ones(1,1191);
NG=ones(1,1191);                                                           % Ngain=1e-20;% Nloss=1e-20;% NG=1;
sigma=0.27;                                                                %���ϵ��,��ֵΪ0.2ʱ�޷��ҵ��㶨�����,���ֹ���0.5,sigmaΪ0.5,NsΪ10ʱcircleΪ100�����ҵ�ƽ��ֵ0.0854
t=sqrt(1-sigma^2);         %����ϵ��

for circle=1:200
%         P1 = (abs(t*sqrt(P1)+i*sigma*sqrt(P)))^2;
%         Ac =t*Ac;
%         As =t*As+i*sigma*As0;
for cnt = 1:Ns
        I1 = P1/Aeff;                                                      %����ǿ�� %N1=1e-6*(1/(1-exp(-T/tau)))*B*PW/(2*En)*I1.^2;%1e-6*(tau/T + 1/(2+T/tau))*B*PW/(2*En)*I1.^2;
        N1 = 1e-6*B*tau*I1.^2/(2*En);                                      %�����������ܶ�
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for NF
        P2=P1;                                                             %����ǿ��
            As2=As;
            Ac2=Ac; 
         NetGains2=ones(1,1191);   
        for dnt=cnt:Ns
            I2 = P2/Aeff;                                                  %���ֹ�ǿ��     %N2=1e-6*(1/(1-exp(-T/tau)))*B*PW/(2*En)*I2.^2;%1e-6*(tau/T + 1/(2+T/tau))*B*PW/(2*En)*I2.^2;
            N2 = 1e-6*B*tau*I2.^2/(2*En);                                  %�����������ܶ�
            dk2= dB+2*ga*P2;                                               %��λʧ����
            g2 = sqrt((ga*P2)^2-(dk2./2).^2);                              %����ϵ��
            Asin2=As2;                                                     %�źŹ����
            Acin2=Ac2;                                                     %��Ƶ�����
            cAcin2=conj(Acin2);                                            %������Ƶ�����
            %%%%%%%%%%%%%%%%%%%%%���㽻����λ���ƺ��Ĳ���Ƶ���start
            PhiP2=exp(-i*ga*P2*dz);
            a2 = cosh(g2.*dz) + i*0.5*(dk2./g2).*sinh(dz*g2);
            b2 = i*((ga*P2)./g2).*exp(i*2*ga*P2*dz).*sinh(dz*g2);
             ca2 = conj(a2);
             cb2 = conj(b2);
             As2 = PhiP2.*(a2.*Asin2 + b2.*cAcin2);
             cAc2 = PhiP2.*(cb2.*Asin2 + ca2.*cAcin2);
             Ac2 = conj(cAc2);
             Gain3=(abs(As2)./abs(Asin2)).^2;
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
             loss2=1.45e-17*N2;                                            %���ֹ������������������
             loss21=1.45e-17*N2*(Ls./Lp).^2;                               %�źŹ������������������
             loss22=2.45e-17*N2*(Lc./Lp).^2;                               %��Ƶ�������������������
             LinLoss3=exp(-dz/2*(100*loss2 + 1*att));                      %�����ʽ�б��ֹ��������������պ����������
             LinLoss31 = exp(-dz/2*(100*loss21 + 1*att));                  %�����ʽ���źŹ��������������պ����������
             LinLoss32 = exp(-dz/2*(100*loss22 + 1*att));                  %�����ʽ����Ƶ���������������պ����������
             LinLoss4 = exp(-dz*(100*loss2 + att));                        %���Ժ����������������
             NonlinLoss2 = exp(-2*dz*B*I2);                                %˫�������������
             NetGains2 = NetGains2.*Gain3.*LinLoss4.*NonlinLoss2;          %������
             P2  = P2*(LinLoss4 * 1./(1+I2*B*dz));
             Ac2 = (LinLoss32*exp(-1/2*2*dz*B*I2)).*Ac2;                   %��������
             As2 = (LinLoss31*exp(-1/2*2/2*dz*B*I2)).*As2;                 %��������
        end      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%���㽻����λ���ƺ��Ĳ���Ƶ���start
        dk   = dB + 2*ga*P1;
        g    = sqrt((ga*P1)^2-(dk./2).^2);
        Asin = As;
        Acin = Ac;
        cAcin= conj(Ac);
        PhiP = exp(-i*ga*P1*dz);
        a    = cosh(g.*dz) + i*0.5*(dk./g).*sinh(dz*g);
        b    = i*((ga*P1)./g).*exp(i*2*ga*P1*dz).*sinh(dz*g);
        ca   = conj(a);
        cb   = conj(b);
        As   = PhiP.*(a.*Asin + b.*cAcin);
        cAc  = PhiP.*(cb.*Asin + ca.*cAcin);
        Ac   = conj(cAc);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end
        loss  = 1.45e-17*N1;                                               %���ֹ������������������
        loss1 = 1.45e-17*N1*(Ls./Lp).^2;                                   %�źŹ������������������
        loss2 = 1.45e-17*N1*(Lc./Lp).^2;                                   %��Ƶ�������������������
        LinLosss = exp(-dz/2*(100*loss1 + att));                           %�����ʽ���źŹ��������������պ����������
        LinLossc = exp(-dz/2*(100*loss2 + att));                           %�����ʽ����Ƶ���������������պ����������
        LinLoss = exp(-dz/2*(100*loss + att));                             %�����ʽ�б��ֹ��������������պ����������
        LinLoss2 = exp(-dz*(100*loss + att));                              %���Ժ����������������
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for NF
        Gain1 = (abs(As)./abs(Asin)).^2;
        NonlinLoss = exp(-2*dz*B.*I1); %�����������
        Ngain = Ngain+NetGains2.*log(Gain1);
        Nloss = Nloss+NetGains2.*(-log(LinLoss2.*NonlinLoss));
        NG = NG.*Gain1.*LinLoss2.*NonlinLoss;       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        P1 = P1*(LinLoss2 * 1./(1+I1*B*dz));
        Ac = (LinLossc*exp(-1/2*2*dz*B*I1)).*Ac;
        As = (LinLosss*exp(-1/2*2*dz*B*I1)).*As;
 
end
         As = exp(i*2*pi*neff*L./Ls).*As;
         Ac = exp(i*2*pi*neff*L./Lc).*Ac;
         P1 = (abs(exp(i*2*pi*neff*L./Lp).*sqrt(P1))).^2;
         
         Asout = t*As0+i*sigma*As;
         Gain5 =(abs(Asout)./abs(As0)).^2;
         LNetG1 = 10*log10(Gain5);
         Gain6 =(abs(As)./abs(As0)).^2;
         LNetG6 = 10*log10(Gain6);
         
         Acout=i*sigma*Ac;
         Gainc   = (abs(Acout)./abs(As0)).^2;
         LNetGc=10*log10(Gainc);
         Gainc6   = (abs(Ac)./abs(As0)).^2;
         LNetGc6=10*log10(Gainc6);
%          plot(1e6*Lc,LNetGc);
         
        P1 = (abs(t*sqrt(P1)+i*sigma*sqrt(P)))^2;
        Ac =t*Ac;
        As =t*As+i*sigma*As0;
        x=angle(As);
        LNF=10*log10((Gain5+Ngain+Nloss)./Gain5);
        LNFc=10*log10((Gainc+Ngain+Nloss)./Gainc);
        plot(1e6*Lc,LNetGc);
%         plot(1e6*Ls,LNetG1,1e6*Lc,LNFc,1e6*Lc,LNF);
%         Result = [1e6*Ls; LNetG1]';
%         save (FName, 'Result', '-ascii');
end





