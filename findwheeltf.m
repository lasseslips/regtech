data = readtable('wheellog.txt');
data = fillmissing(data,'nearest');

t = table2array(data(:,1));     % time stamps
u_L = table2array(data(:,8));   % left motor voltage
u_R = table2array(data(:,9));   % right motor voltage
v_L = table2array(data(:,10));  % left motor velocity 
v_R = table2array(data(:,11));  % right motor velocity
T_s = t(10) - t(9);              % compute sampling time


idx_out = find(u_L >= 3,1,'first');

v_L = v_L(idx_out:end);
v_R = v_R(idx_out:end);
u_L = u_L(idx_out:end);
u_R = u_R(idx_out:end);
u_L = u_L - u_L(1);
v_L = v_L - v_L(1);


t = t(idx_out:end);
t = t - t(1);


vel = 1/2*(v_L + v_R);          % average wheel velocity
volt = 1/2*(u_L + u_R);         % average voltage

idd_L = iddata(v_L,u_L,T_s);    % Prepare left wheel data
idd_R = iddata(v_R,u_R,T_s);
idd = iddata(vel,volt,T_s);
Gvel_tf = tfest(idd,2,1)

G_wu_L = tfest(idd_L,2,0);      % Identify left wheel transfer function (2 poles, no zeros - you can try other combinations)
G_wu_R = tfest(idd_R,2,0); 


kp_R = 15;
kp_L = 15;
w_B_R = 28;
w_B_L = 2;

G_filter_R = tf([w_B_R],[1,w_B_R]);
G_filter_L = tf([w_B_L],[1,w_B_L]);


G_ol_R = G_wu_R*kp_R*G_filter_R;
G_ol_L = G_wu_L*kp_L*G_filter_L;
G_cl_R = (kp_R*G_wu_R)/(1+G_ol_R);
G_cl_L = (kp_L*G_wu_L)/(1+G_ol_L);





