%code to generate LHS and PRCC

clear all
clc
runs=1000;


%PRCC_var={'rho','rq24','w','p','time_v'};% use if waning is included

PRCC_var={'rho','rq24','p','time_v'};% use if waning is not included


par1 = xlsread('par1.xls');%data from fitting

%generate LHS Matrix

rho_LHS=LHS_Call(0.0001,0.5, 1, 0 ,runs,'unif');
rq24_LHS=LHS_Call(0.0001,0.5, 1, 0 ,runs,'unif');
% w_LHS=LHS_Call(0.00001,1/(0.5*365), 1/(0.25*365),0 ,runs,'unif'); % %uncomment if waning is included
p_vacc_LHS=LHS_Call(0.1,0.5, 0.9, 0 ,runs,'unif');
p_time_LHS=LHS_Call(21,100,200, 0 ,runs,'unif');

%  LHSmatrix=[rho_LHS rq24_LHS w_LHS p_vacc_LHS  p_time_LHS]; %uncomment if waning is included
LHSmatrix=[rho_LHS rq24_LHS p_vacc_LHS  p_time_LHS]; %uncomment if waning is not included

options_ode = odeset('NonNegative',1:10,'RelTol', 1e-8, 'AbsTol', 1e-8);

%initial values
Im0=par1(28);
H0=37;
D0=0;
R0=4;
E0=par1(1);
Ia0=par1(2);
Iq0=par1(3);
S0=	2956024-E0-Ia0-Im0-Iq0-R0-H0-D0;
Im0_cumu=150;
start_IV=[S0, E0, Ia0,Im0,Iq0,H0,D0,R0, Im0_cumu,D0];


for x=1:runs
    f=@ODE_before_vaccine;
    f2=@ ODE_after_vaccine;
    paramset=LHSmatrix(x,:);
    
    %uncomment if waning is included
    
    %      rho=paramset(1);
    % rq24=paramset(2);
    % w=paramset(3);
    % p_vacc=paramset(4);
    % time_v=paramset(5);
    
    %uncomment if waning is not included
    rho=paramset(1);
    rq24=paramset(2);
    p_vacc=paramset(3);
    time_v=paramset(4);
    
    
    par4=0.7; %efficacy
    
    [T1,Y1] = ode23(@(t,y) f(t,y,paramset,x,runs,par1,par4),0:1/24:273,start_IV,options_ode);
    start_IV2=Y1(end,:);
    par6=Y1(end,1);
    [T2,Y2] = ode23(@(t,y) f2(t,y,paramset,x,runs,par1,par4,par6),273:1/24:806,start_IV2,options_ode);
    
    
    
    
    z1= Y2(end,9); %cumulative cases
    z2= Y1(end,7); %cumulative deaths
    
    %save solutions for each run
    cumu_cases(:,x)=(z1);
    cumu_deaths(:,x)=(z2);
end
y_var={'Cumulative Cases'};
PRCC_PLOT2(LHSmatrix,cumu_cases,PRCC_var,y_var)



y_var={'Cumulative Deaths'};
PRCC_PLOT2(LHSmatrix,cumu_cases,PRCC_var,y_var)


