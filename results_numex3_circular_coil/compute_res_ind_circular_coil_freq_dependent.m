function [inductance] = compute_res_ind_circular_coil_freq_dependent(rad_wire,rad_loop,freq)

cr_const=2*rad_wire/(2*rad_loop);

% see the condition below for this formula working correctly
if (cr_const > 0.2)
    %disp('formula is not gonna work!!!')
    %disp(['cr constant:::',cr_const])
    %disp('It should be smaller than 0.2')
else
    %disp('formula is gonna work!!!')
end


rad_wire_cm=rad_wire/1e-2;
rad_loop_cm=rad_loop/1e-2;


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

%sig_xx=sig_xx;
inductance=0.01257*rad_loop_cm*...
    (2.303*log10(16*rad_loop_cm/(2*rad_wire_cm))-2+1*sig_xx);

inductance=inductance*1e-6;