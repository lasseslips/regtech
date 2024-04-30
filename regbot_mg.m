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
RA = 3.3/2;    % ohm (2 motors)
JA = 1.3e-6*2; % motor inertia
LA = 6.6e-3/2; % rotor inductor (2 motors)
BA = 3e-6*2;   % rotor friction
Kemf = 0.0105; % motor constant
Km = Kemf;
% køretøj
NG = 9.69; % gear
WR = 0.03; % wheel radius
Bw = 0.155; % wheel distance
% 
% model parts used in Simulink
mmotor = 0.193;   % total mass of motor and gear [kg]
mframe = 0.32;    % total mass of frame and base print [kg]
mtopextra = 0.97 - mframe - mmotor; % extra mass on top (charger and battery) [kg]
mpdist =  0.10;   % distance to lit [m]
% disturbance position (Z)
pushDist = 0.1; % relative to motor axle [m]

%% wheel velocity controller (no balance) PI-regulator
% sample (usable) controller values
Kpwv = 15;     % Kp
tiwv = 0.05;   % Tau_i
Kffwv = 0;     % feed forward constant
startAngle = 10;  % tilt in degrees at time zero
twvlp = 0.005;    % velocity noise low pass filter time constant (recommended)

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
%% Bodeplot
%h = figure(100)
%bode(Gwv)
%grid on
%title('Transfer function from motor voltage to velocity')


%%tilt

iostilt(1) = linio(strcat(model, '/Transfer Fcn'),1,'openinput');
iostilt(2) = linio(strcat(model, '/tilt_angle'),1,'openoutput');
setlinio(model,iostilt);
systilt = linearize(model,iostilt,op);
[numtilt, dentilt] = ss2tf(systilt.A,systilt.B,systilt.C,systilt.D);
Gtilt = minreal(tf(numtilt,dentilt));
%post integrator
kppost = -1.56;
tipost = 0.1; 
Gtilt_post = Gtilt*kppost*tf([tipost,1],[tipost,0]);

%tilt controller
w = logspace(-2,2,3000);

alpha = 0.1;
Ni = 3;
phasemargin = 60;
[wc, Kp_tilt, taui,taud,ok] = findpid(Gtilt_post,phasemargin,Ni,alpha,w);
Gd_tilt = tf([taud 1], [alpha*taud 1]);
Gi_tilt = tf([taui 1], [taui 0]);
G_ol_tilt = Kp_tilt*Gd_tilt*Gi_tilt*Gtilt_post;
G_controller_tilt = Kp_tilt*Gd_tilt*Gi_tilt;
G_cl_tilt = (G_ol_tilt)/(G_ol_tilt+1);


%%velocity controller
Gvel = tf([5736],[1 313.5 1.305e04])*G_cl_tilt;

alpha = 5;
Ni = 3;
phasemargin = 60;
[wc, Kp_vel, taui,taud,ok] = findpid(Gvel,phasemargin,Ni,alpha,w);
Gd_vel = tf([taud 1], [alpha*taud 1]);
Gi_vel = tf([taui 1], [taui 0]);
G_ol_vel = Kp_vel*Gd_vel*Gi_vel*Gvel;
G_controller_vel = Kp_vel*Gd_vel*Gi_vel;
G_cl_vel = (G_ol_vel)/(G_ol_vel+1);


%%position controller
Gpos = G_cl_vel * tf([1],[1 0]);

alpha = 0.01;
Ni = 5;
phasemargin = 70;
[wc, Kp_pos, taui,taud,ok] = findpid(Gpos,phasemargin,Ni,alpha,w);
Gd_pos = tf([taud 1], [alpha*taud 1]);
Gi_pos = tf([taui 1], [taui 0]);
G_ol_pos = Kp_pos*Gd_pos*Gi_pos*Gpos;
G_controller_pos = Kp_pos*Gd_pos*Gi_pos;
G_cl_pos = (G_ol_pos)/(G_ol_pos+1);


