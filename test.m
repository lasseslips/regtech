data = readtable('test.txt');
data = fillmissing(data,'nearest');

t = table2array(data(:,1));     % time stamps
T_s = t(10) - t(9);              % compute sampling time

tilt = table2array(data(:,18));
cntrlbal = table2array(data(:,20));


plot(t,cntrlbal);
