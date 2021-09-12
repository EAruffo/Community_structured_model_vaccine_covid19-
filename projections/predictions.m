%%%this code is to run predictions


clear all
clc


par1 = xlsread('par1.xls');%parameters from fitting
par2_vec=[0.1, 0.3, 0.6, 0.9]; % vacciantion coverages
par3_vec=[60; 180]; %time reaching the peak of the vaccination function
par3_lab=[60;180]; %labels for the peak of the vaccination function
par4_vec=[0.7, 0.9]; %efficacy of vaccine
par4_lab=[70, 90];% labels for efficacy of vaccine
par5_vec=[0,1/(0.25*365), 1/(0.5*365), 1/(365) ]; % waning rates
par5_lab=[0,90, 180, 365 ]; %labels for waning rates


for wan=1:length(par5_vec)
    
    par5= par5_vec(wan);
    
    for k=1:length(par4_vec)
        par4=par4_vec(k);
        
        for j=1:length(par3_vec)
            par3=par3_vec(j);
            
            for i=1:length(par2_vec)
                par2=par2_vec(i);
                
                
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
                
                %the running time is divided into two phases: before and
                %after vaccination
                
                [T1,Y1] = ode23(@ODE_before_vaccination,0:1/24:273,start_IV,options_ode,par1,par2,par3,par4,par5);
                start_IV2=Y1(end,:);
                par6=Y1(end,1); % susceptibles when the vaccine starts
                [T2,Y2] = ode23(@ODE_after_vaccination,273:1/24:806,start_IV2,options_ode,par1,par2,par3,par4,par5,par6);
                
                
                %plot of cumulative deaths
                
                figure(j)
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                plot(0:273,Y1(1:24:end,7),'*r')
                hold on
                plot(273:806,Y2(1:24:end,7),'LineWidth',3)
                ax = gca; % current axes
                ax.FontSize = 24;
                ax.FontWeight = 'bold';
                title(([' peak=' num2str(par3_lab(j)) '  ']),'FontSize',24,'FontWeight','Bold');
                xlabel('Time ')
                ylabel('Cumulative deaths')
                xlim([0 807])
                %change the title if using different phases
                %change the xtiks depending if the reopening occurs in May
                %or December
                %xticks([0 273 624 806 ])
                %xticklabels({'Mar17,20',' Dec15,20 ', 'Dec 1,21',' Jun 1,22'})
                xticks([0 273 410 806 ])
                xticklabels({'Mar17,20',' Dec15,20 ', 'May 1,21',' Jun 1,22'})
                
                %plot of cumulative cases
                
                figure(j+1)
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                plot(0:273,Y1(1:24:end,9),'*r')
                hold on
                plot(273:806,Y2(1:24:end,9),'LineWidth',3)
                ax = gca; % current axes
                ax.FontSize = 24;
                ax.FontWeight = 'bold';
                title(([' peak=' num2str(par3_lab(j)) '  ']),'FontSize',24,'FontWeight','Bold');
                xlabel('Time ')
                ylabel('Cumulative cases')
                xlim([0 807])
                title((['Phase 2: reopening on May 1, waning=' num2str(par5_lab(wan)) '  ']),'FontSize',24,'FontWeight','Bold');%change the title if using different phases

                %change the xtiks depending if the reopening occurs in May
                %or December
                %xticks([0 273 624 806 ])
                %xticklabels({'Mar17,20',' Dec15,20 ', 'Dec 1,21',' Jun 1,22'})
                xticks([0 273 410 806 ])
                xticklabels({'Mar17,20',' Dec15,20 ', 'May 1,21',' Jun 1,22'})
                
                title((['Phase 2: reopening on May 1, waning=' num2str(par5_lab(wan)) '  ']),'FontSize',24,'FontWeight','Bold');%change the title if using different phases
                
            end
            close all
            
            
        end
        
        
    end
end