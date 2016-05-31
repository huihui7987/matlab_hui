%
%高斯脉冲的微分
%函数2:"cp0702_Gaussian_derivatives"
%先将高斯脉冲及其前15阶导函数的幅度归一化，
%然后将它们同时画出。
%该函数仅有一输值'alpha'（即脉冲形成因α），
%没有输出值。
%执行函数2的命令如下：
%cp0702_Gaussian_derivatives(alpha);
%
% FUNCTION 2 : "cp0702_Gaussian_derivatives"
%
% Analysis of waveforms of the Gaussian pulse and its first
% 15 derivatives
%
% The pulse amplitude is set to 'A'
% 'smp' samples of the Gaussian pulse are considered in
% the time interval 'Tmax - Tmin'
%
% The function receives as input the value of the shape
% factor 'alpha'
%
% The function plots in a 4 x 4 grid the waveform of the
% Gaussian pulse and of its first 15 derivatives for the
% 'alpha' received as input
% 
% Programmed by Luca De Nardis
 
%function cp0702_Gaussian_derivatives(alpha)
 
% -----------------------------------------------
% Step Zero - Input parameters and initialization
% -----------------------------------------------
alpha=0.714e-9;
A = 1;                            % pulse amplitude [V]
smp = 1024;                       % number of samples
Tmin = -4e-9;                     % lower time limit
Tmax = 4e-9;                      % upper time limit
 
t=linspace(Tmin,Tmax,smp);        % initialization of the
                                  % time axis
pulse=-A*exp(-2*pi*(t/alpha).^2); % pulse waveform
                                  % definition
 
F=figure(1);
set(F,'Position',[100 190 850 450]);
subplot(4,4,1);
PT=plot(t,pulse);
axis([-2e-9 2e-9 -1 1]);
set(gca,'XTick',0);
set(gca,'XTickLabel',{
 
});
 
for i=1:15
    % determination of the i-th derivative
    derivative(i,:) = ...
       cp0702_analytical_waveforms(t,i,alpha);
 
    % amplitude normalization of the i-th derivative
    derivative(i,:) = derivative(i,:) / ...
       max(abs(derivative(i,:)));
 
% -----------------------------------------------
% Step One - Graphical output
% -----------------------------------------------
 
    subplot(4,4,i+1);
    PT=plot(t,derivative(i,:));
    axis([-2e-9 2e-9 -1 1]);
    if(i < 12)
        set(gca,'XTick',0);
        set(gca,'XTickLabel',{
 
});
    end
    if(mod(i,4) ~= 0)
        set(gca,'YTickLabel',{
 
});
    end
end
h = axes('Position',[0 0 1 1],'Visible','off');
set(gcf,'CurrentAxes',h);
text(.5,0.02,'Time[s]','FontSize',12)
text(0.05,0.45,'Amplitude [V]','FontSize',12,...
   'Rotation', 90);