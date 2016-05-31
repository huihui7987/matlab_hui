function re = gausplus(a,b)
%程序输入参数a抽样点数、b脉冲类型
    samples=a;  type=b;  t=((1:samples)'-1)/samples-0.5;
    %定义高斯脉宽(波形约为半宽)Tp=pi*tao且fc=1/Tp,若一阶导、二阶导数波形
    %的中心频率为fc1和fc2,则有fc1=fc/sqrt(2);fc2=fc;
    Tp=0.1;   tao=Tp/pi;  
    switch(type)
    %高斯脉冲与其傅氏变换形式相同(钟形函数)
    case 0, gaus=exp(-t.^2/tao.^2); re=gaus/norm(gaus); 
    %高斯单脉冲与其傅氏变换形式相同,但谱幅度函数是对时域函数取绝对值,变成偶函数     
    %脉宽为Tp的高斯脉冲的一阶导数ygmono对应的中心频率为fc1=fc/sqrt(2)
    %按道理来讲tao取相等值时,一阶导数也应与函数库中的gmonopuls(t,fc1)相对应
    %但却与函数库中的gmonopuls(t,fc)相对应,理解成二者的表达式是不完全同的
    %gmonopuls(..)并不是高斯波形直接求导的结果,经过了中心频率归一化
    case 1, gausmono = -2*t./tao^2.*exp(-t.^2/tao.^2); re=gausmono/norm(gausmono); 
    %gausmono2=gmonopuls(t,fc);%=-2*t./tao^2.*exp(1/2*(-t.^2/tao.^2)); 
    %高斯二阶导数波形
    case 2,  gausdoub=(-2+4.*t.^2./tao^2)./tao^2.*exp(-t.^2/tao.^2); re=gausdoub/norm(gausdoub); 
    %高斯波形五阶导数
    case 5,  sigma=0.08448;   gausfive=-(-15.*t./sigma+10.*(t./sigma).^3-(t./sigma).^5).*exp(-(t./sigma).^2./2); re=gausfive/norm(gausfive);
    %Scholtz单脉冲，其中心频率为fp=sqrt(2/pi)/tao,比如tao=1/pi有fp=5
otherwise, scholtz=sqrt(8/3/tao)*(1-4*pi*(t/tao).^2).*exp(-2*pi*(t/tao).^2); re=scholtz/norm(scholtz);
    end  