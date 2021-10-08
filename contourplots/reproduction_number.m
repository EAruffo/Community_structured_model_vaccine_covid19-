%%this codes generates contourplots for the reproduction number
clear all;
close all;
clc;

N=200; % points on the grid

%parameters
alpha=(1/4);
b=0.8;
gamma_ia=0.07;
gamma_mr=(1/14);
ch=2.4;
gammaq24=1; 
gammaq48=1/2; 
gamma_mh=1/6;

%parameters to infectigate
p_value=[0,0.3,0.9];
rho=[0:1/200:1];
r24=[0:1/200:1];

%uncomment the phase to use for the reproduction number
%Phase 1
p2=0.144;
p1=1-p2;
cc=4.378645821;
cw=6.839067896;
betah_f=6.72E-08;
bhl=betah_f-betah_f*0.1;
bhu=betah_f+betah_f*0.1;
step_h=(bhu-bhl)/N;
betah=[bhl:step_h:bhu];


%Phase 2
% p2=0.0663;
% p1=1-p2;
% cc=7.421244069;
% cw=8.914671221;
% betah_f=2.65e-08;
%  bhl=betah_f-betah_f*0.1;
%  bhu=betah_f+betah_f*0.1;
%  step_h=(bhu-bhl)/N;
%  betah=[bhl:step_h:bhu];

% %Phase 3
% p2=0.0469;
% p1=1-p2;
% cc=7.905791217;
% cw=9.937479748;
% betah_f=6.86e-08;
%  bhl=betah_f-betah_f*0.1;
%  bhu=betah_f+betah_f*0.1;
%  step=(bhu-bhl)/N;
%  betah=[bhl:step:bhu];

%phase 4
% p2=0.0794;
% p1=1-p2;
% cc=6.423286525;
% cw=9.709616148;
% betah_f=5.8e-08;
%  bhl=betah_f-betah_f*0.1;
%  bhu=betah_f+betah_f*0.1;
%  step=(bhu-bhl)/N;
%  betah=[bhl:step:bhu];

G=1;% constant to stay at home

[X,Y,Z] = meshgrid(rho,r24,betah);


for j=1:1%length(p_value)
    
    
    S0=2954489-2954489*p_value(j)*0.7; %proportion of susceptibles after vaccination 70% effective
    
R0_A= (1-b)/gamma_ia;
R0_I=G*b./((1-X).*(p1*gamma_mr+p2*gamma_mh) + X.*Y*gammaq24 + X.*gammaq48.*(1-Y));

R=Z.*ch*S0.*(R0_A+R0_I);


%90% effective
S02=2954489-2954489*p_value(j)*0.9; %proportion of susceptibles after vaccination 70% effective
R0_A2= (1-b)/gamma_ia;
R0_I2=G*b./((1-X).*(p1*gamma_mr+p2*gamma_mh) + X.*Y*gammaq24 + X.*gammaq48.*(1-Y));

R2=Z.*ch*S02.*(R0_A2+R0_I2);

cvals = linspace(0,2,N);

Sx = [];
Sy = [];
Sz = [betah(1),betah(25),betah(50),betah(100),betah(N)];



figure(1)
subplot(2,3,j)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

csh=contourslice(X,Y,Z,R,Sx,Sy,Sz,cvals);
set(csh,'linewidth',3)
colb=colorbar;
zlim([bhl max(bhu)])
 Az = -17;
 El = 30;
view(Az, El);
grid on

xlabel('\rho','FontSize',24,'FontWeight','Bold');
ylabel('r_{q24}','FontSize',24,'FontWeight','Bold');
zlabel('\beta_h','FontSize',24,'FontWeight','Bold');
title((['p = ' num2str(p_value(j)) '']),'FontSize',14,'FontWeight','Bold')
set(gca,'FontSize',24)



hold on
subplot(2,3,j+3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

csh2=contourslice(X,Y,Z,R2,Sx,Sy,Sz,cvals);
set(csh2,'linewidth',3)
%caxis([0 2])
colb=colorbar;
set(gca,'FontSize',24)
zlim([bhl max(bhu)])
 Az = -17;
 El = 30;
view(Az, El);
grid on
xlabel('\rho','FontSize',24,'FontWeight','Bold');
ylabel('r_{q24}','FontSize',24,'FontWeight','Bold');
zlabel('\beta_h','FontSize',24,'FontWeight','Bold');
title((['p = ' num2str(p_value(j)) '']),'FontSize',24,'FontWeight','Bold')

a=axes;
tt=title('R_0 in the household in Phase 1, efficacy 70%-90%','FontSize',13,'FontWeight','Bold')
set(a,'Visible', 'off')
            set(tt,'Visible', 'on')
            set(tt,'Position',get(tt,'Position')+[0 .04 0])


end
