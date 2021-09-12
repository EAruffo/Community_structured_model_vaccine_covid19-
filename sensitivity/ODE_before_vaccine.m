
function dy= ODE_before_vaccine(t,y,paramset,x,runs,par1,par4)

% max_vacc=paramset(4); %uncomment if waning is included
max_vacc=paramset(3); %uncomment if waning is not included

S0=0;
% time_P=paramset(5); %uncomment if waning is included
time_P=paramset(4); %uncomment if waning is not included

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
f_v=[0.6601 0.6985 0.2261 0.7075];


hours=t-floor(t);
if t>= 0 & t<76
    p=0;
    w=0;
    rho=0.635;
    y1_t=0*ones(size(t));
    Y_v=y1_t;
    
    mu_H=par1(9);
    gamma_qh=par1(8);
    rq=rq24I;
    p2=p2_v(1);
    p1=p1_v(1);
    f=f_v(1);
    q=q_v(1);

    if  hours>= 9/24  & hours < 15/24
        G=0;
        beta=  par1(6);
        c=par1(5); 
    elseif hours >= 15/24  & hours < 17/24
        G=0;
        beta= par1(6);
        c=par1(4);
        
    else
        G=1;
        beta= par1(7);
        c=c_h; 
        
    end
    
elseif t>=76 & t<160
    p=0;
    w=0;
    y1_t=0*ones(size(t));
    Y_v=y1_t;
    rho=0.635;
    
    mu_H=par1(15);
    gamma_qh=par1(14);
    rq=rq24II;
    p2=p2_v(2);
    p1=p1_v(2);
    f=f_v(2);
    q=q_v(2);
    if hours>= 9/24  & hours < 15/24
        G=0;
        beta= par1(12);
        c=par1(11); %c_w
    elseif hours >= 15/24 & hours < 19/24
        G=0;
        beta= par1(12);
        c=par1(10); %c_c
    else
        G=1;
        beta= par1(13);
        c=c_h; %c_h
        
    end
    
    
elseif t>=160 & t<238
    p=0;
    w=0;
    y1_t=0*ones(size(t));
    Y_v=y1_t;
    rho=0.635;
    
    p2=p2_v(3);
    p1=p1_v(3);
    f=f_v(3);
    q=q_v(3);
    
    mu_H=par1(21);
    gamma_qh=par1(20);
    rq=rq24III;
    
    if hours >= 9/24 & hours < 15/24
        G=0;
        beta= par1(18);
        c=par1(17); 
    elseif hours >= 15/24 & hours < 17/24
        G=0;
        beta= par1(18);
        c=par1(16); 
    else
        G=1;
        beta=  par1(19);
        c=c_h; 
        
    end

elseif t>=238 & t<261
    p=0;
    w=0;
    y1_t=0*ones(size(t));
    Y_v=y1_t;
    rho=0.635;
    
    p2=p2_v(4);
    p1=p1_v(4);
    f=f_v(4);
    q=q_v(4);
    mu_H=par1(27);
    gamma_qh=par1(26);
    
    rq=rq24III;
    
    if hours >= 9/24 & hours < 15/24
        G=0;
        beta= par1(24);
        c=par1(23); 
    elseif hours >= 15/24 & hours < 17/24
        G=0;
        beta= par1(24);
        c=par1(22);
    else
        G=1;
        beta=  par1(25);
        c=c_h; 
        
    end

else
    p=0;
    w=0;
    y1_t=0*ones(size(t));
    Y_v=y1_t;
   rho=paramset(1);
    
%     %PHASE1
%  mu_H=par1(9);
%         gamma_qh=par1(8);
%         p2=p2_v(1);
%         p1=p1_v(1);
%         f=f_v(1);
%         q=q_v(1);
%     rq=paramset(2);
%     
%         if  hours>= 9/24  & hours < 15/24
%             G=0;
%             beta=  par1(6);
%             c=par1(5); 
%         elseif hours >= 15/24  & hours < 17/24
%             G=0;
%             beta= par1(6);
%             c=par1(4);
%     
%         else
%             G=1;
%             beta= par1(7);
%             c=c_h; 
%     
%         end
        

%PHASE 2
mu_H=par1(15);
         gamma_qh=par1(14);
rq=paramset(2);
          p2=p2_v(2);
        p1=p1_v(2);
        f=f_v(2);
        q=q_v(2);
        
        
        if hours>= 9/24  & hours < 15/24
            G=0;
            beta= par1(12);
            c=par1(11); 
        elseif hours >= 15/24 & hours < 19/24
            G=0;
            beta= par1(12);
            c=par1(10); 
        else
            G=1;
            beta= par1(13);
            c=c_h; 
    
         end

end


S=y(1); E=y(2); Ia=y(3); Im=y(4);  Iq=y(5);  H=y(6); D= y(7);  R= y(8);

dy=[-p*eff*Y_v*S0-beta.*c.*S.*(Ia+G.*Im) + w*R;%S
    beta.*c.*S.*(Ia+G.*Im)-alpha.*E; %E
    (1-b).*alpha.*E-gamma_ia.*Ia; %Ia
    b*alpha*E-(1-rho)*(p1*gamma_mr + p2*gamma_mh)*Im - ...
    rho*rq*gammaq24*Im- rho*(1-rq)*gammaq48*Im; %I
    rho*rq*gammaq24*Im + rho*(1-rq)*gammaq48*Im - (1-q)*gamma_qh*Iq-q*gamma_mr*Iq; %Iq
    (1-rho).*p2.*gamma_mh.*Im + (1-q).*(gamma_qh).*Iq-(1-f)*gamma_H.*H- f*mu_H*H; %H
    f*mu_H*H; %D
    p*eff*Y_v*S0+(1-rho).*p1.*gamma_mr.*Im + q.*(gamma_mr).*Iq+ (1-f).*gamma_H.*H+gamma_ia.*Ia- w*R; %R
    b*alpha*E; f*mu_H*H]; %cumulative c, d


end

