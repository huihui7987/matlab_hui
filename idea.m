clear all
lamda=1562.7e-9:1e-12:1563.1e-9;
v=(3e8./lamda)-(3e8./1562.9e-9);

%%%%%%%%%
%阶数
%%%%%%%%%
n1=0.5;
n2=1;
n3=1.5;
%p1=exp(-1i*pi/2);
%p2=exp(1i*pi/2);
%t=2*pi*v;

%%%%%%%%%%理想透射谱%%%%%%%%%%%%
H1=mapminmax((abs(1i*2*pi*v)).^(2*n1),0,1);
H2=mapminmax((abs(1i*2*pi*v)).^(2*n2),0,1);
H3=mapminmax((abs(1i*2*pi*v)).^(2*n3),0,1);
plot(v,H1,'r'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,H2,'g'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
plot(v,H3,'b'); xlabel('Frequency(Hz）');ylabel('Intensity Transmission');hold on;
legend('0.5','1','1.5')
figure;
%%%%%%%%理想相位响应%%%%%%%%%%%%%%
p1=(n1*pi/2).*(v>0)+(-n1*pi/2).*(v<0);
p2=(n2*pi/2).*(v>0)+(-n2*pi/2).*(v<0);
p3=(n3*pi/2).*(v>0)+(-n3*pi/2).*(v<0);

plot(v,p1,'r','linewidth',2.5);hold on;
plot(v,p2,'g','linewidth',2.5); hold on;
plot(v,p3,'b','linewidth',2.5);hold on;
% plot(v,p1,'r');hold on;
% plot(v,p2,'g'); hold on;
% plot(v,p3,'b');hold on;

figure;
%%%%%%%%%%%%%%%%理想传输函数%%%%%%%%%%%%%%%%%

Hv1=((1i*2*pi*v)).^(n1);
Hv2=((1i*2*pi*v)).^(n2);
Hv3=((1i*2*pi*v)).^(n3);
subplot(1,3,1)
plot(v,Hv1,'r'); xlabel('Frequency(Hz）');ylabel('yyyyy');
subplot(1,3,2)
plot(v,Hv2,'r'); xlabel('Frequency(Hz）');ylabel('yyyyy');
subplot(1,3,3)
plot(v,Hv3,'r'); xlabel('Frequency(Hz）');ylabel('yyyyy');
figure;

H1=(abs(2*pi*v)).*(n1)*exp(1i*n1*pi/2).*(v>0)+(abs(2*pi*v)).*(n1)*exp(-1i*n1*pi/2).*(v<0);
H2=(abs(2*pi*v)).*(n2)*exp(1i*n2*pi/2).*(v>0)+(abs(2*pi*v)).*(n2)*exp(-1i*n2*pi/2).*(v<0);
H3=(abs(2*pi*v)).*(n3)*exp(1i*n3*pi/2).*(v>0)+(abs(2*pi*v)).*(n3)*exp(-1i*n3*pi/2).*(v<0);
subplot(1,3,1)
plot(v,H1,'r'); xlabel('Frequency(Hz）');ylabel('xxx');hold on;
subplot(1,3,2)
plot(v,H2,'g'); xlabel('Frequency(Hz）');ylabel('xxx');hold on;
subplot(1,3,3)
plot(v,H3,'b'); xlabel('Frequency(Hz）');ylabel('xxx');hold on;
legend('0.5','1','1.5')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%高斯函数经过微分
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%高斯函数%%%%%%%












