

r=0.83;
tao=0.801;
R=50E-6;
lamda=1.438:0.001:1.440;
neff=3.179;

L = 2*pi*R;
phi = mod(L*neff./lamda*2*pi,2*pi);%~~~~-
Ta = exp(1i*(pi+phi)).*(tao-r.*exp(-1i*phi))./(1-r.*tao.*exp(1i*phi));
T = (abs(Ta)).^2;%~~~~-
%T = (tao^2-2*r*tao*cos(phi)+r^2)./(1-2*r*tao*cos(phi)+r^2*tao^2);
PHI = angle(Ta);%~~~~-
% if r<=tao
%     PHI = PHI+(PHI<0)*2*pi ;
% end;

plot(lamda,PHI)





