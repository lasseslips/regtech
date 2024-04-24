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
Kpwv = 20;     % Kp
tiwv = 0.01;   % Tau_i
Kffwv = 0;     % feed forward constant
startAngle = 0;  % tilt in degrees at time zero
twvlp = 0.005;    % velocity noise low pass filter time constant (recommended)

%% Estimate transfer function for base system using LINEARIZE
% Motor volatge to wheel velocity (wv)
load_system(model);
open_system(model);
% define points in model
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
[mag phase] = bode(Gtilt_post,w);
alpha = 0.1;
Ni = 3;
phasemargin = 60;
phi_d = rad2deg(asin((1-alpha)/(1+alpha)));
phi_i = rad2deg(atan2(-1,Ni));
pc = phasemargin - 180 - phi_d - phi_i;
n = find(phase > pc, 1, 'last');
wc = w(n);
tau_d = 1/(sqrt(alpha)*wc);
tau_i = Ni/wc;
Gd_tilt = tf([tau_d 1],[alpha*tau_d 1]);
Gi_tilt = tf([tau_i 1], [tau_i 0]);
[magc, phasec] = bode(Gtilt_post*Gi_tilt*Gd_tilt,wc);
Kp_tilt = 1/magc;
Gtilt_ol = Gtilt_post*Gi_tilt*Gd_tilt*Kp_tilt;
Gtilt_cl = minreal((Gtilt_ol)/(Gtilt_ol+1));
Gtilt_controller = Gi_tilt*Gd_tilt*Kp_tilt;
