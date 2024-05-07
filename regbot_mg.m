%% Simscape multibody model og Regbot in balance
% initial setup with motor velocity controller 
% this is intended as simulation base for balance control.
%
close all
clear

%% Simulink model name
model='regbot_1yayamg';

%% parameters for REGBOT
% motor
%given values
% RA = 3.3/2;    % ohm (2 motors)
% JA = 1.3e-6*2; % motor inertia
% LA = 6.6e-3/2; % rotor inductor (2 motors)
% BA = 3e-6*2;   % rotor friction
% Kemf = 0.0105; % motor constant
% Km = Kemf;

%measured values
% RA = (6+12)/4;
% JA = (0.00093+0.00065)/(4*100);
% LA = 6.6e-3/2;
% BA = (0.00000383+0.00000058);
% Kemf = (0.00812+0.00805)/2;
% Km = Kemf;

%test
RA = (5.5)/2;
JA = (0.00093+0.00065)/(4*1000);
LA = 6.6e-3/2;
BA = (0.00000383+0.00000058);
Kemf = (0.00812+0.00805)/2;
Km = Kemf;

% køretøj
NG = 9.69; % gear
WR = 0.03; % wheel radius
Bw = 0.155; % wheel distance



% model parts used in Simulink
mmotor = 0.193;   % total mass of motor and gear [kg]
mframe = 0.32;    % total mass of frame and base print [kg]
mtopextra = 0.97 - mframe - mmotor; % extra mass on top (charger and battery) [kg]
mpdist =  0.10;   % distance to lit [m]
% disturbance position (Z)
pushDist = 0.1; % relative to motor axle [m]

%% wheel velocity controller (no balance) PI-regulator
% sample (usable) controller values
Kpwv = 10.13;     % Kp
tiwv = 0.2727;   % Tau_i
Kffwv = 0;     % feed forward constant
startAngle = 40;  % tilt in degrees at time zero
twvlp = 0.005;    % velocity noise low pass filter time constant (recommended)

%initial values for simulink
kppost = -1;
taui_post = 1;
taud_tilt = 1;
G_controller_velu_num = 1;
G_controller_velu_den = 1;
G_controller_tilt_num = 1;
G_controller_tilt_den = 1;
taufilter = 1;
G_controller_vel_num = 1;
G_controller_vel_den = 1;
G_controller_pos_num = 1;
G_controller_pos_den = 1;

%% Estimate transfer function for base system using LINEARIZE
% Motor volatge to wheel velocity (wv)
load_system(model);
open_system(model);
% define points inG model
ios(1) = linio(strcat(model,'/Limit9v'),1,'openinput');
ios(2) = linio(strcat(model, '/wheel_vel_filter'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gwv = minreal(tf(num, den))

%Gwv = tf([1.1775 1.68],[1 7.334 7.141]);

w = logspace(-2,2,3000);

alpha = 5;
Ni = 4;
phasemargin = 80;

[wc, Kp_velu, taui,taud,ok] = findpid(Gwv,phasemargin,Ni,alpha,w);
Gd_velu = tf([taud 1], [alpha*taud 1]);
Gi_velu = tf([taui 1], [taui 0]);

G_ol_velu = minreal(Kp_velu*Gd_velu*Gi_velu*Gwv);
%Kp_velu = 10;
G_controller_velu = Kp_velu*Gi_velu;
G_controller_velu_num = G_controller_velu.Numerator{1};
G_controller_velu_den = G_controller_velu.Denominator{1};
G_cl_velu = minreal((G_ol_velu)/(G_ol_velu+1));





Gol = Gwv*Kpwv*tf([tiwv 1],[tiwv 0]);
Gcl = (Gol)/(Gol + 1);



%%tilt

iostilt(1) = linio(strcat(model, '/veluRef'),1,'openinput');
iostilt(2) = linio(strcat(model, '/tilt_angle'),1,'openoutput');
setlinio(model,iostilt);
systilt = linearize(model,iostilt,op);
[numtilt, dentilt] = ss2tf(systilt.A,systilt.B,systilt.C,systilt.D);
Gtilt = minreal(tf(numtilt,dentilt));
%%post integrator
kppost = -1;
[gain,freq] = getPeakGain(-Gtilt);
taui_post = 1/freq;
Gi_post = tf([taui_post 1], [taui_post 0]);
Gtilt_post = Gtilt*kppost*Gi_post;
Gtilt_post_cl = minreal(Gtilt_post/(Gtilt_post + 1));

%tilt controller
alpha = 0.01;
Ni = 4;
phasemargin = 80;

[wc, Kp_tilt, taui_tilt,taud_tilt,ok] = findpid(Gtilt_post_cl,phasemargin,Ni,alpha,w);
Gd_tilt = tf([taud_tilt 1], [alpha*taud_tilt 1]);
Gi_tilt = tf([taui_tilt 1], [taui_tilt 0]);

G_ol_tilt = minreal(Kp_tilt*Gd_tilt*Gi_tilt*Gtilt_post_cl);
G_controller_tilt = Kp_tilt*Gi_tilt*Gd_tilt;
G_controller_tilt_num = G_controller_tilt.Numerator{1};
G_controller_tilt_den = G_controller_tilt.Denominator{1};
G_cl_tilt = minreal((G_ol_tilt)/(G_ol_tilt+1));

%%tilt prefilter
[mag phase wout] = bode(G_cl_tilt);
[gain freq] = findpeaks(squeeze(mag));
taufilter = 1/wout(freq(end));
%taufilter = 1;
G_tilt_filter = tf([1],[taufilter 1]);
G_cl_tilt_filter = G_tilt_filter*G_cl_tilt;


%%velocity controller
ios(1) = linio(strcat(model,'/tiltRef'),1,'openinput');
ios(2) = linio(strcat(model, '/wheel_vel_filter'),1,'openoutput');
setlinio(model,ios);
op = [0];
sys = linearize(model,ios,op);
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gvel = minreal(tf(num, den));


alpha = 0.9;
Ni = 5;
phasemargin = 70;

[wc, Kp_vel, taui,taud,ok] = findpid(Gvel,phasemargin,Ni,alpha,w);
Gd_vel = tf([taud 1], [alpha*taud 1]);
Gi_vel = tf([taui 1], [taui 0]);

G_ol_vel = Kp_vel*Gd_vel*Gi_vel*Gvel;
G_controller_vel = Kp_vel*Gd_vel*Gi_vel;
G_controller_vel_num = G_controller_vel.Numerator{1};
G_controller_vel_den = G_controller_vel.Denominator{1};
G_cl_vel = (G_ol_vel)/(G_ol_vel+1);


%%position controller
ios(1) = linio(strcat(model,'/velRef'),1,'openinput');
ios(2) = linio(strcat(model, '/xposRef'),1,'openoutput');
setlinio(model,ios);
op = [0];
sys = linearize(model,ios,op);
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gpos = minreal(tf(num, den));






alpha = 1;
Ni = 15;
phasemargin = 70;
[wc, Kp_pos, taui,taud,ok] = findpid(Gpos,phasemargin,Ni,alpha,w);
Gd_pos = tf([taud 1], [alpha*taud 1]);
Gi_pos = tf([taui 1], [taui 0]);
G_ol_pos = Kp_pos*Gd_pos*Gi_pos*Gpos;
G_controller_pos = Kp_pos*Gd_pos*Gi_pos;
G_controller_pos_num = G_controller_pos.Numerator{1};
G_controller_pos_den = G_controller_pos.Denominator{1};
G_cl_pos = (G_ol_pos)/(G_ol_pos+1);


Kp_velu
Gd_velu
Gi_velu

Gi_post
Kp_tilt
Gd_tilt
Gi_tilt
G_tilt_filter

Kp_vel
Gd_vel
Gi_vel

Kp_pos
Gd_pos
Gi_pos