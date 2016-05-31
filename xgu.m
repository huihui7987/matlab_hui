function re = gausplus(a,b)
%�����������a����������b��������
    samples=a;  type=b;  t=((1:samples)'-1)/samples-0.5;
    %�����˹����(����ԼΪ���)Tp=pi*tao��fc=1/Tp,��һ�׵������׵�������
    %������Ƶ��Ϊfc1��fc2,����fc1=fc/sqrt(2);fc2=fc;
    Tp=0.1;   tao=Tp/pi;  
    switch(type)
    %��˹�������丵�ϱ任��ʽ��ͬ(���κ���)
    case 0, gaus=exp(-t.^2/tao.^2); re=gaus/norm(gaus); 
    %��˹���������丵�ϱ任��ʽ��ͬ,���׷��Ⱥ����Ƕ�ʱ����ȡ����ֵ,���ż����     
    %����ΪTp�ĸ�˹�����һ�׵���ygmono��Ӧ������Ƶ��Ϊfc1=fc/sqrt(2)
    %����������taoȡ���ֵʱ,һ�׵���ҲӦ�뺯�����е�gmonopuls(t,fc1)���Ӧ
    %��ȴ�뺯�����е�gmonopuls(t,fc)���Ӧ,���ɶ��ߵı��ʽ�ǲ���ȫͬ��
    %gmonopuls(..)�����Ǹ�˹����ֱ���󵼵Ľ��,����������Ƶ�ʹ�һ��
    case 1, gausmono = -2*t./tao^2.*exp(-t.^2/tao.^2); re=gausmono/norm(gausmono); 
    %gausmono2=gmonopuls(t,fc);%=-2*t./tao^2.*exp(1/2*(-t.^2/tao.^2)); 
    %��˹���׵�������
    case 2,  gausdoub=(-2+4.*t.^2./tao^2)./tao^2.*exp(-t.^2/tao.^2); re=gausdoub/norm(gausdoub); 
    %��˹������׵���
    case 5,  sigma=0.08448;   gausfive=-(-15.*t./sigma+10.*(t./sigma).^3-(t./sigma).^5).*exp(-(t./sigma).^2./2); re=gausfive/norm(gausfive);
    %Scholtz�����壬������Ƶ��Ϊfp=sqrt(2/pi)/tao,����tao=1/pi��fp=5
otherwise, scholtz=sqrt(8/3/tao)*(1-4*pi*(t/tao).^2).*exp(-2*pi*(t/tao).^2); re=scholtz/norm(scholtz);
    end  