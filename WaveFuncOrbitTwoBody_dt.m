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


%% use wavefunction to split the parallel and pependicular motion
% calculate the orbit over time.
% Caution only x-y plane motion is simulated, z has to be zero. also only apogee at
% x-axis is simulated. Random orientation and starting point have not been
% coded.

clear;
% close all;

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

% caculate the total energy and a factor for both objects
miu1 = m1^3/(m1+m2)^2;
miu2 = m2^3/(m1+m2)^2;
r10 = sqrt((x10-xc)^2 + (y10-yc)^2 + (z10-zc)^2);
r20 = sqrt((x20-xc)^2 + (y20-yc)^2 + (z20-zc)^2);
E1 = -G*m1*miu2/r10 + 0.5*m1*(vx10^2+vy10^2+vz10^2);
E2 = -G*m2*miu1/r20 + 0.5*m2*(vx20^2+vy20^2+vz20^2);
RE1 = -G*m1*miu2/2/E1;
RE2 = -G*m2*miu1/2/E2;

vcr1 = cross([x10-xc, y10-yc, z10-zc], [vx10, vy10, vz10]);
vcr1 = sqrt(dot(vcr1, vcr1)); % vr*r1
a10 = vcr1^2/G/miu2/RE1;
vcr2 = cross([x20-xc, y20-yc, z20-zc], [vx20, vy20, vz20]);
vcr2 = sqrt(dot(vcr2, vcr2));
a20 = vcr2^2/G/miu1/RE2;


vr10 = vy10; % same should be perp to radius. (real perpendicular to radius)
vi10 = 0; %(imaginary part value parallel to radius)
vr20 = -vy20; 
vi20 = 0;

%% use Newton's law to calculate locations and speed
%initial condition
data1 = zeros(length(t), 6);  % t, x, y, z, vr, vi 
data1(:,1) = t; 
data1(1, 2) = x10;
data1(1, 3) = y10;
data1(1, 4) = z10;
data1(1, 5) = vr10;
data1(1, 6) = vi10;


data2 = zeros(length(t), 6);  % t, x, y, z, vr, vi
data2(:,1) = t; 
data2(1, 2) = x20;
data2(1, 3) = y20;
data2(1, 4) = z20;
data2(1, 5) = vr20;
data2(1, 6) = vi20;


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


    
    % calculate m1 orbit
    x11 = data1(i-1, 2); y11 = data1(i-1, 3); z11 = data1(i-1, 4);
    vr11 = data1(i-1, 5); vi11 = data1(i-1, 6);
    r11 = sqrt((x11-xc)^2 + (y11-yc)^2 + (z11-zc)^2);
    theta11 = atan2d(y11-yc, x11-xc)*pi()/180;
    thetad1 = vr11*dt/r11;
    theta12 = theta11 + thetad1;
    if y11 < 0
        approach1 = true;
    else
        approach1 = false;
    end
    if approach1
        r12 = r11*cos(thetad1) - vi11*dt; % - OR + ???
    else
        r12 = r11/cos(thetad1) + vi11*dt; % - OR + ???
    end
    x12 = r12*cos(theta12) + xc;
    y12 = r12*sin(theta12) + yc;
    z12 = 0; % no z motion for now
    

    vr12 = abs(sqrt(a10*G*miu2*RE1)/r12); %some case near 0 has a leak and cause v to be imaginary, so add abs.
    vi12 = abs(sqrt(-G*miu2/RE1 + 2*G*miu2/r12 - a10*G*miu2*RE1/r12^2));
    
    data1(i, 2) = x12; % m1 state
    data1(i, 3) = y12;
    data1(i, 4) = z12;
    data1(i, 5) = vr12;
    data1(i, 6) = vi12;
    
    
    % calculate m2 orbit
    x21 = data2(i-1, 2); y21 = data2(i-1, 3); z21 = data2(i-1, 4);
    vr21 = data2(i-1, 5); vi21 = data2(i-1, 6); 
    r21 = sqrt((x21-xc)^2 + (y21-yc)^2 + (z21-zc)^2);
    
    theta21 = atan2d(y21-yc, x21-xc)*pi()/180;
    thetad2 = vr21*dt/r21;
    theta22 = theta21 + thetad2;
    if y21 >= 0
        approach2 = true;
    else
        approach2 = false;
    end
    if approach2
        r22 = r21*cos(thetad2) - vi21*dt; % - OR + ???
    else
        r22 = r21/cos(thetad2) + vi21*dt; % - OR + ???
    end
    x22 = r22*cos(theta22) + xc;
    y22 = r22*sin(theta22) + yc;
    z22 = 0; % no z motion for now
    
    vr22 = abs(sqrt(a20*G*miu1*RE2)/r22); %some case near 0 has a leak and cause v to be imaginary, so add abs.
    vi22 = abs(sqrt(-G*miu1/RE2 + 2*G*miu1/r22 - a20*G*miu1*RE2/r22^2));
    
    data2(i, 2) = x22; % m1 state
    data2(i, 3) = y22;
    data2(i, 4) = z22;
    data2(i, 5) = vr22;
    data2(i, 6) = vi22; 
    
          
    
    dataC(i, 2) = x12 + (x22-x12)*m2/(m1+m2); %CoM, should not change but check for error.
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


for i = 1: length(t)
    plot(data1(i,2), data1(i,3), 'ro');
    plot(data2(i,2), data2(i,3), 'bo');
    plot(dataC(i,2), dataC(i,3), 'k+', 'markersize', 10, 'linewidth', 2); %CoM % should not change but check for error.
%     plot(dataC(i,5), dataC(i,6), 'kx'); %CoE 
    pause(0.001);
end
hold off

