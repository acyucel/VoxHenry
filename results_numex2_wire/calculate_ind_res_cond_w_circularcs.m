function [resistance_dc, resistance_ac, inductance_ac] = calculate_ind_res_cond_w_circularcs(freq,len_wire,rad_wire,se)

%clc; close all; clear all;
%freq=3e0; len_wire=50.0e-6; rad_wire=5.0e-6; se=5.8e7; % conductivity of interconnect

mu0=4*pi*1e-7;
rhoa_res=1/se; % resistivity of wire

A_area=pi*rad_wire^2;

% DC resistance
resistance_dc=rhoa_res*len_wire/A_area;

% AC resistance

skin_depth=sqrt(rhoa_res/(pi*freq*mu0));

% calculation of effective area with various methods
% (http://chemandy.com/calculators/round-wire-ac-resistance-calculator.htm)

% method 1
% A_area_eff=skin_depth*pi*(2*rad_wire);
% remark on method 1: This method makes the used cross sectional area 
% too large from high frequencies down to the point 
% where skin depth becomes about half of the radius of the conductor, 
% at which point the inaccuracies increase and eventually the calculated 
% used area becomes bigger than the actual conductor. 
% Making the calculation method only an approximation and then 
% only usable when r >> ?


% method 2
% A_area_eff=(pi*rad_wire^2)-(pi*(rad_wire-skin_depth)^2);

% remark on method 2: This method is more accurate the first method 
% when r >> ? but becomes highly inaccurate below the point 
% where d/? = ? and can have huge positive or negative swings in value.

% method 3

skin_depth_prime=skin_depth*(1-exp(-rad_wire/skin_depth));

zz_dum=0.62006*rad_wire/skin_depth;
yy_dum=0.189774/(1+0.272481*(zz_dum^1.82938-zz_dum^-0.99457)^2)^1.0941;
A_area_eff=pi*(2*rad_wire*skin_depth_prime-skin_depth_prime^2)*(1+yy_dum);

resistance_ac=rhoa_res*len_wire/A_area_eff;

resistance_dc = rhoa_res*len_wire/A_area;

% the formula in 1937 report
%dia_wire=2*rad_wire
%xx_a=0.1405*2*rad_wire_cm*sqrt(mu0*freq*5.8e11);
rad_wire_cm=rad_wire/1e-2;
len_wire_cm=len_wire/1e-2;
dia_wire_cm=2*rad_wire_cm;

xx_a=0.1071*2*rad_wire_cm*sqrt(freq)*(1.68/1.72);

fl_lin_interp=1;

if (fl_lin_interp == 1)
    %linear_interp(y0,y1,x0,x1,xx);
    if(xx_a>=0 && xx_a<=0.5)
        sig_xx=linear_interp(0.250,0.250,0,0.5,xx_a);
    elseif (xx_a > 0.5 && xx_a<=1.0)
        sig_xx=linear_interp(0.250,0.249,0.5,1.0,xx_a);
    elseif (xx_a > 1.0 && xx_a<=1.5)
        sig_xx=linear_interp(0.249,0.247,1.0,1.5,xx_a);
    elseif (xx_a > 1.5 && xx_a<=2.0)
        sig_xx=linear_interp(0.247,0.240,1.5,2.0,xx_a);
    elseif (xx_a > 2.0 && xx_a<=2.5)
        sig_xx=linear_interp(0.240,0.228,2.0,2.5,xx_a);
    elseif (xx_a > 2.5 && xx_a<=3.0)
        sig_xx=linear_interp(0.228,0.211,2.5,3.0,xx_a);
    elseif (xx_a > 3.0 && xx_a<=3.5)
        sig_xx=linear_interp(0.211,0.191,3.0,3.5,xx_a);
    elseif (xx_a > 3.5 && xx_a<=4.0)
        sig_xx=linear_interp(0.191,0.1715,3.5,4.0,xx_a);
    elseif (xx_a > 4.0 && xx_a<=4.5)
        sig_xx=linear_interp(0.1715,0.154,4.0,4.5,xx_a);
    elseif (xx_a > 4.5 && xx_a<=5.0)
        sig_xx=linear_interp(0.154,0.139,4.5,5.0,xx_a);
    elseif (xx_a > 5.0 && xx_a<=6.0)
        sig_xx=linear_interp(0.139,0.116,5.0,6.0,xx_a);
    elseif (xx_a > 6.0 && xx_a<=7.0)
        sig_xx=linear_interp(0.116,0.100,6.0,7.0,xx_a);
    elseif (xx_a > 7.0 && xx_a<=8.0)
        sig_xx=linear_interp(0.100,0.088,7.0,8.0,xx_a);
    elseif (xx_a > 8.0 && xx_a<=9.0)
        sig_xx=linear_interp(0.088,0.078,8.0,9.0,xx_a);
    elseif (xx_a > 9.0 && xx_a<=10.0)
        sig_xx=linear_interp(0.078,0.070,9.0,10.0,xx_a);
    elseif (xx_a > 10.0 && xx_a<=12.0)
        sig_xx=linear_interp(0.070,0.059,10.0,12.0,xx_a);
    elseif (xx_a > 12.0 && xx_a<=14.0)
        sig_xx=linear_interp(0.059,0.050,12.0,14.0,xx_a);
    elseif (xx_a > 14.0 && xx_a<=16.0)
        sig_xx=linear_interp(0.050,0.044,14.0,16.0,xx_a);
    elseif (xx_a > 16.0 && xx_a<=18.0)
        sig_xx=linear_interp(0.044,0.039,16.0,18.0,xx_a);
    elseif (xx_a > 18.0 && xx_a<=20.0)
        sig_xx=linear_interp(0.039,0.035,18.0,20.0,xx_a);
    elseif (xx_a > 20.0 && xx_a<=25.0)
        sig_xx=linear_interp(0.035,0.028,20.0,25.0,xx_a);
    elseif (xx_a > 25.0 && xx_a<=30.0)
        sig_xx=linear_interp(0.028,0.024,25.0,30.0,xx_a);
    elseif (xx_a > 30.0 && xx_a<=40.0)
        sig_xx=linear_interp(0.024,0.0175,30.0,40.0,xx_a);
    elseif (xx_a > 40.0 && xx_a<=50.0)
        sig_xx=linear_interp(0.0175,0.014,40.0,50.0,xx_a);
    elseif (xx_a > 50.0 && xx_a<=60.0)
        sig_xx=linear_interp(0.014,0.012,50.0,60.0,xx_a);
    elseif (xx_a > 60.0 && xx_a<=70.0)
        sig_xx=linear_interp(0.012,0.010,60.0,70.0,xx_a);
    elseif (xx_a > 70.0 && xx_a<=80.0)
        sig_xx=linear_interp(0.010,0.009,70.0,80.0,xx_a);
    elseif (xx_a > 80.0 && xx_a<=90.0)
        sig_xx=linear_interp(0.009,0.008,80.0,90.0,xx_a);
    elseif (xx_a > 90.0 && xx_a<=100.0)
        sig_xx=linear_interp(0.008,0.007,90.0,100.0,xx_a);
    else
        sig_xx=0.0;
    end
    
elseif (fl_lin_interp == 0)
    
    if(xx_a>=0 && xx_a<=0.5)
        sig_xx=0.250;
    elseif (xx_a > 0.5 && xx_a<=1.0)
        sig_xx=0.249;
    elseif (xx_a > 1.0 && xx_a<=1.5)
        sig_xx=0.247;
    elseif (xx_a > 1.5 && xx_a<=2.0)
        sig_xx=0.240;
    elseif (xx_a > 2.0 && xx_a<=2.5)
        sig_xx=0.228;
    elseif (xx_a > 2.5 && xx_a<=3.0)
        sig_xx=0.211;
    elseif (xx_a > 3.0 && xx_a<=3.5)
        sig_xx=0.191;
    elseif (xx_a > 3.5 && xx_a<=4.0)
        sig_xx=0.1715;
    elseif (xx_a > 4.0 && xx_a<=4.5)
        sig_xx=0.154;
    elseif (xx_a > 4.5 && xx_a<=5.0)
        sig_xx=0.139;
    elseif (xx_a > 5.0 && xx_a<=6.0)
        sig_xx=0.116;
    elseif (xx_a > 6.0 && xx_a<=7.0)
        sig_xx=0.100;
    elseif (xx_a > 7.0 && xx_a<=8.0)
        sig_xx=0.088;
    elseif (xx_a > 8.0 && xx_a<=9.0)
        sig_xx=0.078;
    elseif (xx_a > 9.0 && xx_a<=10.0)
        sig_xx=0.070;
    elseif (xx_a > 10.0 && xx_a<=12.0)
        sig_xx=0.059;
    elseif (xx_a > 12.0 && xx_a<=14.0)
        sig_xx=0.050;
    elseif (xx_a > 14.0 && xx_a<=16.0)
        sig_xx=0.044;
    elseif (xx_a > 16.0 && xx_a<=18.0)
        sig_xx=0.039;
    elseif (xx_a > 18.0 && xx_a<=20.0)
        sig_xx=0.035;
    elseif (xx_a > 20.0 && xx_a<=25.0)
        sig_xx=0.028;
    elseif (xx_a > 25.0 && xx_a<=30.0)
        sig_xx=0.024;
    elseif (xx_a > 30.0 && xx_a<=40.0)
        sig_xx=0.0175;
    elseif (xx_a > 40.0 && xx_a<=50.0)
        sig_xx=0.014;
    elseif (xx_a > 50.0 && xx_a<=60.0)
        sig_xx=0.012;
    elseif (xx_a > 60.0 && xx_a<=70.0)
        sig_xx=0.01;
    elseif (xx_a > 70.0 && xx_a<=80.0)
        sig_xx=0.009;
    elseif (xx_a > 80.0 && xx_a<=90.0)
        sig_xx=0.008;
    elseif (xx_a > 90.0 && xx_a<=100.0)
        sig_xx=0.007;
    else
        sig_xx=0.0;
    end
    
end


% uncomment the following
inductance_ac=0.002*len_wire_cm*(2.303*log10(4*len_wire_cm/dia_wire_cm)-1+1*sig_xx+(dia_wire_cm/(2*len_wire_cm)));
inductance_ac=(inductance_ac)*1e-6; 

% % compute inductance
% % formula in http://chemandy.com/calculators/round-wire-inductance-calculator.htm
% % convert the quantities from m to cm
% rad_wire_cm=rad_wire/1e-2;
% len_wire_cm=len_wire/1e-2;
% dia_wire_cm=2*rad_wire_cm;
% 
% xx=2*pi*rad_wire_cm*sqrt(2.0*mu0*freq/(se)) % simens per cm
% 
% T_xx=sqrt((0.873011+0.00186128*xx)/(1.0-0.278381*xx+0.127964*(xx^2)));
% 
% inductance_ac_uH=0.002*len_wire_cm*(log(4.0*len_wire_cm/dia_wire_cm) ...
%     -1.0 + (dia_wire_cm/(2.0*len_wire_cm)) + (1.0/4.0 *  T_xx));
% 
% 
% inductance_ac=inductance_ac_uH*1e-6;

% inductance_ac=0.002*len_wire_cm*(2.303*log10(4*len_wire_cm/dia_wire_cm)-1+1*sig_xx+(dia_wire_cm/(2*len_wire_cm)));
% inductance_ac=(inductance_ac)*1e-6; %+0.22e-11


% new inductance calculation with the formula 1 in 
% https://en.wikipedia.org/wiki/Inductor#Inductance_formulas
% rosa's formula
% const1 = mu0/(2*pi);
% 
% const2 =(len_wire)*log((1/rad_wire)*(len_wire+sqrt(len_wire^2+rad_wire^2)));
% 
% const3 = -sqrt(len_wire^2+rad_wire^2)+rad_wire;
% 
% const4 = len_wire / (4+rad_wire*sqrt(2*se*2*pi*freq*mu0));
% 
% inductance_ac = 0;
% 
% inductance_ac=const1*(const2+const3+const4);


