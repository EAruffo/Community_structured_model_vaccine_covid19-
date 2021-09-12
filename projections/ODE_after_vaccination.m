
function dy= ODE_after_vaccination(t,y,par1,par2,par3,par4,par5,par6)


max_vacc=par2; %this is for vaccine coverage proportion p

S_T0=par6;
%here we define the shape to inform vaccination
time_P=par3; %time at which peak is reached
max_value=1/time_P;


%parameters
eff=par4;
alpha=(1/4);
b=0.8;
gamma_ia=0.07;
gamma_mr=(1/14);
c_h=2.4;
gammaq24=1;
gammaq48=1/2;
rq24I=0.5;
rq24II=0.45;
rq24III=0.37;
gamma_mh=1/6;
gamma_H=1/10;
p2_v=[0.144 0.0663 0.0469 0.0794];
p1_v=1-p2_v;
q_v=p1_v;
rho=0.635;
f_v=[0.6601 0.6985 0.2261 0.7075];

hours=t-floor(t); %define the daily hours

%define the shape of the vacciantion function 

if t>=273 & t<273+time_P
    
    y2_t=(t)*max_value/time_P - max_value*273/time_P;
elseif t>=273+time_P & t<273+2*time_P
    y2_t=-(t)*max_value/time_P + max_value*(273+2*time_P)/time_P;
else
    y2_t=0*ones(size(t));
    
end

Y_v=y2_t;


%simulations

if t>=273 & t<624 %NOTE: change the econd value, reopening in December 624 or May 410
    
    
    p=max_vacc;
    w=par5;
    
    %    %PHASE1
%      mu_H=par1(9);
%         gamma_qh=par1(8);
%         rq=rq24I;
%         p2=p2_v(1);
%         p1=p1_v(1);
%         f=f_v(1);
%         q=q_v(1);

%         if  hours>= 9/24  & hours < 15/24 %working time
%             G=0;
%             beta=  par1(6);
%             c=par1(5); 
%         elseif hours >= 15/24  & hours < 17/24 % community time
%             G=0;
%             beta= par1(6);
%             c=par1(4);
%     
%         else %home time
%             G=1;
%             beta= par1(7);
%             c=c_h; 
%     
%         end
    
    %PHASE2
%     mu_H=par1(15);
%          gamma_qh=par1(14);
%         rq=rq24II;
%           p2=p2_v(2);
%         p1=p1_v(2);
%         f=f_v(2);
%         q=q_v(2);

%         if hours>= 9/24  & hours < 15/24 %working time
%             G=0;
%             beta= par1(12);
%             c=par1(11);
%         elseif hours >= 15/24 & hours < 19/24 %community time
%             G=0;  
%             beta= par1(12);
%             c=par1(10); 
%         else %home time
%             G=1;
%             beta= par1(13);
%             c=c_h; 
%     
%          end
    
    %PHASE3
%     p2=p2_v(3);
%     p1=p1_v(3);
%     f=f_v(3);
%     q=q_v(3);
%     mu_H=par1(21);
%     gamma_qh=par1(20);
%     rq=rq24III;
%     
%     if hours >= 9/24 & hours < 15/24 %working time
%         G=0;
%         beta= par1(18);
%         c=par1(17); 
%     elseif hours >= 15/24 & hours < 17/24 %community time
%         G=0;
%         beta= par1(18);
%         c=par1(16); 
%     else %home time
%         G=1;
%         beta=  par1(19);
%         c=c_h;
%         
%     end
    %PHASE4
     p2=p2_v(4);
        p1=p1_v(4);
        f=f_v(4);
        q=q_v(4);
        mu_H=par1(27);
         gamma_qh=par1(26);
    
        rq=rq24III;
    
        if hours >= 9/24 & hours < 15/24 %working time
            G=0; 
            beta= par1(24);
            c=par1(23); 
        elseif hours >= 15/24 & hours < 17/24 %community time
            G=0;
            beta= par1(24);
            c=par1(22); 
        else %home time
            G=1;
            beta=  par1(25);
            c=c_h; 
    
        end
   
  %reopening with highest transmission and phase specific parameters  
  %uncomment the phase to simulate
else 
    
    p=max_vacc;
    w=par5;
    
    
%PHASE1
%      mu_H=par1(9);
%         gamma_qh=par1(8);
%         rq=rq24I;
%         p2=p2_v(1);
%         p1=p1_v(1);
%         f=f_v(1);

    %PHASE2
%       mu_H=par1(15);
%          gamma_qh=par1(14);
%         rq=rq24II;
%           p2=p2_v(2);
%         p1=p1_v(2);
%         f=f_v(2);
%         q=q_v(2);

    %PHASE3
%     p2=p2_v(3);
%     p1=p1_v(3);
%     f=f_v(3);
%     q=q_v(3);
%     mu_H=par1(21);
%     gamma_qh=par1(20);
%     rq=rq24III;

    %PHASE4
         p2=p2_v(4);
        p1=p1_v(4);
        f=f_v(4);
        q=q_v(4);
        mu_H=par1(27);
         gamma_qh=par1(26);
        rq=rq24III;
    
    
    
    if  hours>= 9/24  & hours < 15/24 %working time
        G=0;
        beta=  par1(6);
        c=par1(17);
    elseif hours >= 15/24  & hours < 17/24 %community time
        G=0;
        beta= par1(6);
        c=par1(16);
    else %home time
        G=1;
        beta= par1(19);
        c=c_h;
    end
    
    
    
end

%system

S=y(1); E=y(2); Ia=y(3); Im=y(4);  Iq=y(5);  H=y(6); D= y(7);  R= y(8);

dy=[-p*eff*Y_v*S_T0-beta.*c.*S.*(Ia+G.*Im) + w*R;%S
    beta.*c.*S.*(Ia+G.*Im)-alpha.*E; %E
    (1-b).*alpha.*E-gamma_ia.*Ia; %Ia
    b*alpha*E-(1-rho)*(p1*gamma_mr + p2*gamma_mh)*Im - ...
    rho*rq*gammaq24*Im- rho*(1-rq)*gammaq48*Im; %I
    rho*rq*gammaq24*Im + rho*(1-rq)*gammaq48*Im - (1-q)*gamma_qh*Iq-q*gamma_mr*Iq; %Iq
    (1-rho).*p2.*gamma_mh.*Im + (1-q).*(gamma_qh).*Iq-(1-f)*gamma_H.*H- f*mu_H*H; %H
    f*mu_H*H; %D
    p*eff*Y_v*S_T0+(1-rho).*p1.*gamma_mr.*Im + q.*(gamma_mr).*Iq+ (1-f).*gamma_H.*H+gamma_ia.*Ia- w*R; %R
    b*alpha*E; f*mu_H*H]; %cumulative c, d,






end

