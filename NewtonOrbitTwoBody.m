% copyright 2024 Jixin Chen @ Ohio University, Athens, Ohio
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software. 


%% use Newton's law of gravity F = -GMm/r^2, and F = d(dv/dt)/dt to
% calculate the orbit over time.

clear;
close all;

% initialize parameters
G = 1;
m1 = 1;
m2 = 2;
xc = 0; yc = 0; zc = 0; %center of mass set to be stationary as the inertia reference frame, thus two bodies have correlated freedom.

dt = 0.01;
t = 0:dt:100;
t = t';

%initialize data
x10 = 2;  % initialize state of m1
y10 = 0;
z10 = 0;
vx10 = 0;
vy10 = 0.8;
vz10 = 0;

x20 = (xc-x10)*m1/m2 + xc; % state of m2 is clculatied to keep the center of mass stationary
y20 = (yc-y10)*m1/m2 + yc;
z20 = (zc-z10)*m1/m2 + zc;
vx20 = -vx10*(m1/m2);
vy20 = -vy10*(m1/m2);
vz20 = -vz10*(m1/m2);


r0 = sqrt((x10-x20)^2 + (y10-y20)^2 + (z10-z20)^2);


%% use Newton's law to calculate locations and speed
%initial condition
data1 = zeros(length(t), 7);  % t, x, y, z, vx, vy, vz
data1(:,1) = t; 
data1(1, 2) = x10;
data1(1, 3) = y10;
data1(1, 4) = z10;
data1(1, 5) = vx10;
data1(1, 6) = vy10;
data1(1, 7) = vz10;

data2 = zeros(length(t), 7);  % t, x, y, z, vx, vy, vz
data2(:,1) = t; 
data2(1, 2) = x20;
data2(1, 3) = y20;
data2(1, 4) = z20;
data2(1, 5) = vx20;
data2(1, 6) = vy20;
data2(1, 7) = vz20;

dataC = zeros(length(t), 7); %t, Center of mass (CoM) x, y, z, and Center of Energy (CoE) x, y, z
dataC(:,1) = t; 
dataC(1, 2) = x10 + (x20-x10)*m2/(m1+m2); %CoM
dataC(1, 3) = y10 + (y20-y10)*m2/(m1+m2);
dataC(1, 4) = z10 + (z20-z10)*m2/(m1+m2);
dataC(1, 5) = x10 + (x20-x10)*m1/(m1+m2); %CoE
dataC(1, 6) = y10 + (y20-y10)*m1/(m1+m2);
dataC(1, 7) = z10 + (z20-z10)*m1/(m1+m2);

% calculate the state stepwisely
for i = 2:length(t)
    %read previous data
    x11 = data1(i-1, 2); y11 = data1(i-1, 3); z11 = data1(i-1, 4);
    vx11 = data1(i-1, 5); vy11 = data1(i-1, 6); vz11 = data1(i-1, 7);
    x21 = data2(i-1, 2); y21 = data2(i-1, 3); z21 = data2(i-1, 4);
    vx21 = data2(i-1, 5); vy21 = data2(i-1, 6); vz21 = data2(i-1, 7);
    
    r1 = sqrt((x11-x21)^2 + (y11-y21)^2 + (z11-z21)^2);
    if r1
        Fx1 = G*m1*m2*(x21-x11)/r1^3; Fx2 = G*m1*m2*(x11-x21)/r1^3;
        Fy1 = G*m1*m2*(y21-y11)/r1^3; Fy2 = G*m1*m2*(y11-y21)/r1^3;
        Fz1 = G*m1*m2*(z21-z11)/r1^3; Fz2 = G*m1*m2*(z11-z21)/r1^3;
    else %r1 = 0, singularity
        Fx1 = 0; Fx2 = 0;
        Fy1 = 0; Fy2 = 0;
        Fz1 = 0; Fz2 = 0;
    end
    
    ax1 = Fx1/m1; ay1 = Fy1/m1; az1 = Fz1/m1;
    ax2 = Fx2/m2; ay2 = Fy2/m2; az2 = Fz2/m2;
    
    vx12 = vx11 + ax1*dt; vy12 = vy11 + ay1*dt; vz12 = vz11 + az1*dt;
    vx22 = vx21 + ax2*dt; vy22 = vy21 + ay2*dt; vz22 = vz21 + az2*dt;
    
    x12 = x11 + vx12*dt; y12 = y11 + vy12*dt; z12 = z11 + vz12*dt;
    x22 = x21 + vx22*dt; y22 = y21 + vy22*dt; z22 = z21 + vz22*dt;
        
    
    data1(i, 2) = x12; % m1 state
    data1(i, 3) = y12;
    data1(i, 4) = z12;
    data1(i, 5) = vx12;
    data1(i, 6) = vy12;
    data1(i, 7) = vz12;
    
    data2(i, 2) = x22; % m2 state
    data2(i, 3) = y22;
    data2(i, 4) = z22;
    data2(i, 5) = vx22;
    data2(i, 6) = vy22;
    data2(i, 7) = vz22;
    
    dataC(i, 2) = x12 + (x22-x12)*m2/(m1+m2); %CoM should not change but check for error.
    dataC(i, 3) = y12 + (y22-y12)*m2/(m1+m2);
    dataC(i, 4) = z12 + (z22-z12)*m2/(m1+m2);
    dataC(i, 5) = x12 + (x22-x12)*m1/(m1+m2); %CoE
    dataC(i, 6) = y12 + (y22-y12)*m1/(m1+m2);
    dataC(i, 7) = z12 + (z22-z12)*m1/(m1+m2);

end

figure; hold on
plot3(data1(:,2), data1(:,3), data1(:,4));
plot3(data2(:,2), data2(:,3), data2(:,4));
plot3(dataC(:,2), dataC(:,3), dataC(:,4), 'k+', 'markersize', 10, 'linewidth', 2); %CoM
% plot3(dataC(:,5), dataC(:,6), dataC(:,7),'kx'); %CoE

%% make the movie

figure; hold on
plot(data1(:,2), data1(:,3));
plot(data2(:,2), data2(:,3));
plot(dataC(:,2), dataC(:,3), 'k+', 'markersize', 10, 'linewidth', 2); %CoM

for i = 1:length(t)
    plot(data1(i,2), data1(i,3), 'ro');
    plot(data2(i,2), data2(i,3), 'bo');
    plot(dataC(i,2), dataC(i,3), 'k+', 'markersize', 10, 'linewidth', 2); %CoM
%     plot(dataC(i,5), dataC(i,6), 'kx'); %CoE 
    pause(0.001);
end
hold off
