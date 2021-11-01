%% Tomas Glasberger
% Space vector PWM principle

close all; 
clear all;
%% Basic simulation parameters
t = 0; % actual simulation time
dt = 1e-6; % simulation time step
tmax = 0.04; % maximum simulation time
no_vect = 18; % number of vectors for one period of output signal
 
f_out = 50; % output freuqency
t_pwm = 1/(no_vect*f_out); % corresponding switching frequency
om = 2*pi*f_out; % anfgular velocity
dfi = om*dt;
fi = 0;
t_per = 0*t_pwm; % time in one switching period (in microprocessor corresponds with a counter) 
u_saw=t_pwm; % starting point

fi_per=0;
um = 0.477; % demanded output magnitude
Udc = 1; % dc-link voltage
result = zeros (floor(tmax/dt),9); %simulation results
result_ix = 0; % index in to simulation results array

sector = 0;  % sector of demanded votlage vector position 1-6
t1 = 0;      % time for separate components
t2 = 0;
t0 = 0;


sw_comb = [1,-1,-1; 1,1,-1; -1,1,-1; -1,1,1; -1,-1,1; 1,-1,1]; % switching combination of a converter
V1 = [-1,-1,-1]; % output vectors created using switching combinations
V2 = [-1,-1,-1];
V0 = [-1,-1,-1];
V7 = [1,1,1];
 
u_out = [0, 0, 0]; % inverter phase voltage (ua0, ub0, uc0)
gamma= 0; % angle of the output voltage vector in one sector (0°-60°)
Um=0; % size of the demanded space vector 

uaff = 0; % filtered load phase voltage (corresponds with current waveform for RL load)

ura = 0;
urb=0;
urc=0;

%% Main loop of the simulation
while t<tmax

    ua = um * cos (fi);
    ub = um * cos (fi-2*pi/3);
    uc = um * cos (fi+2*pi/3);
    
    ux = 2./3*(ua - 0.5*ub - 0.5*uc);
    uy = 2./3.*(sqrt(3.)/2*ub - sqrt(3.)/2*uc);
    %% SVPWM
        sector = floor(fi/(pi/3))+1; % find the proper sector where the demanded output voltage vector is placed
        gamma = fi - (sector-1)*pi/3; % angle of the output voltage vector ...in the given sector (0°-60°)
        
        Um = sqrt(ux^2 + uy^2); % size of the demanded vector
          switch sector
             case 1, V1 = sw_comb(1,:);
                     V2 = sw_comb(2,:);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);                                
             case 2, V1 = sw_comb(3,:);
                     V2 = sw_comb(2,:);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);
             case 3, V1 = sw_comb(3,:);
                     V2 = sw_comb(4,:);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);
             case 4, V1 = sw_comb(5,:);
                     V2 = sw_comb(4,:);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);
             case 5, V1 = sw_comb(5,:);
                     V2 = sw_comb(6,:);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);
             case 6, V1 = sw_comb(1,:);
                     V2 = sw_comb(6,:);
                     t2 = abs(Um)/(2/pi*Udc) * t_pwm * sin (pi/3 - gamma)/sin(pi/3);
                     t1 = abs(Um)/(2/pi*Udc) * t_pwm * sin (gamma)/sin(pi/3);
            end;
            t0 = t_pwm - t1 -t2;
        
            %% Calculate corresponding modulation signals
            switch sector
                case 1, ura = t1 +t2 +t0/2;
                        urb = t2+t0/2; 
                        urc = t0/2;
                case 2, ura = t2+t0/2;
                        urb = t1 +t2 + t0/2;
                        urc = t0/2;
                case 3, ura = t0/2;
                        urb = t1 +t2 + t0/2;
                        urc = t2+t0/2;
                case 4, ura = t0/2;
                        urb = t2 + t0/2;
                        urc = t1 + t2 + t0/2;
                case 5, ura = t2 + t0/2;
                        urb = t0/2;
                        urc = t1 + t2 + t0/2;
                case 6, ura = t1 + t2 + t0/2;
                        urb = t0/2;
                        urc = t2+t0/2;
            end;
    if (t_per>=t_pwm) 
            t_per=0; 
 
%             fi_per = fi_per + t_pwm*2*pi*f_out;
%             if fi_per>=2*pi 
%                 fi_per = fi_per - 2*pi; end
        end  
    %simulation of the switching signals and the inverter 
    if (t_per<t0/4) u_out = V0; end
    if ((t_per>=t0/4) && t_per<=(t0/4+t1/2)) u_out = V1; end
    if ((t_per>(t0/4+t1/2)) && (t_per<=(t0/4+t1/2+t2/2))) u_out = V2; end
    if ((t_per>(t0/4+t1/2+t2/2)) && (t_per<=(t0/4+t1/2+t2/2+t0/2))) u_out = V7; end
    if ((t_per>(t0/4+t1/2+t2/2+t0/2)) && (t_per<=(t0/4+t1/2+t2/2+t0/2+t2/2))) u_out = V2; end
    if ((t_per>(t0/4+t1/2+t2/2+t0/2+t2/2)) && (t_per<=(t0/4+t1/2+t2/2+t0/2+t2/2+t1/2))) u_out = V1; end
    if (t_per>(t0/4+t1/2+t2/2+t0/2+t2/2+t1/2)) u_out = V0; end
    
    %% calculation of phase load voltages according to inverter phase voltages
    u_out = u_out*Udc/2;   
    uaf = (2*u_out(1) - u_out(2) - u_out(3))/3;
    ubf = (2*u_out(2) - u_out(1) - u_out(3))/3;
    ucf = (2*u_out(3) - u_out(1) - u_out(2))/3;
  
    uaff = uaff + 0.0005*(uaf - uaff);
    
    uab = uaf - ubf;

    %% Calculation of triangular (saw) signal
    if(t_per==0) 
        u_saw=t_pwm; 
    end  
    if(t_per<t_pwm/2)
        u_saw = u_saw-2*dt;        
    end
    if(t_per>t_pwm/2 && t_per<t_pwm)
        u_saw = u_saw+2*dt;        
    end
    %% Saving results    
    result_ix = result_ix + 1;
    result (result_ix,:) = [t, u_saw, ura,urb,urc, u_out(1), u_out(2), u_out(3),sector];
 
    
    
    %% Go to the next simulation step
    fi = fi + dfi;
    if (fi>=2*pi) 
        fi =fi-2*pi; 
    end
    
    t = t + dt;
    t_per = t_per + dt;
    
    
end
%%
figure;
subplot(5,1,1); plot (result (:,1),result (:,2),result (:,1), result (:,3),...
    result (:,1), result (:,4),...
    result (:,1), result (:,5),...
    'LineWidth',2);grid on;axis([0.00222,0.00338,0,0.0012]);
subplot(5,1,2); plot (result (:,1), result (:,6),'LineWidth',2);grid on;
axis([0.00222,0.00338,-0.51,0.51]);

subplot(5,1,3); plot (result (:,1), result (:,7),'LineWidth',2);grid on;
axis([0.00222,0.00338,-0.51,0.51]);

subplot(5,1,4); plot (result (:,1), result (:,8),'LineWidth',2);grid on;
axis([0.00222,0.00338,-0.51,0.51]);

subplot(5,1,5); plot (result (:,1), result (:,9),'LineWidth',2);grid on;axis([0.00217,0.0033,0,6]);


% figure;
% subplot(5,1,1); plot (result (:,1),result (:,2),result (:,1), result (:,3),...
%     result (:,1), result (:,4),...
%     result (:,1), result (:,5),...
%     'LineWidth',2);grid on;
% subplot(5,1,2); plot (result (:,1), result (:,6),'LineWidth',2);grid on;
% 
% 
% subplot(5,1,3); plot (result (:,1), result (:,7),'LineWidth',2);grid on;
% 
% 
% subplot(5,1,4); plot (result (:,1), result (:,8),'LineWidth',2);grid on;
% 
% 
% subplot(5,1,5); plot (result (:,1), result (:,9),'LineWidth',2);grid on;
% 

%% kresli 1 obrazek
% figure;
% set(gcf,'color',[1,1,1]);
% set(gcf,'Position',[529   759   960   288]);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 32 10]);
% 
% plot (result (:,1), result (:,4)/t_pwm-0.5, 'k','LineWidth',2);grid on; 
% grid on; axis([0.0 0.04 -0.6 0.6]);
% xlabel('x [°]', 'FontSize',16); ylabel ('ms_{a} [-]', 'FontSize',16);
% set(gca, 'FontSize', 16, 'FontName', 'Arial'); print( gcf, '-dpng', '-r300','d:\uz_no_load');

%% FFT
% DESCRIPTIVE TEXT
% usa = (result(:,3));%(result(:,4))/t_pwm-0.5;
% x = fft(usa);
% X = abs(x);
% X = X/(length(usa)/2);
% figure; bar ( X(1:1000));
