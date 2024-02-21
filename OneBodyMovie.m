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
close all
% initialize parameters
GM = 1;
m = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be origin and stationary

%initialize data
dt = 0.001;
t = 0:dt:6;
t = t';
x0 = 1;
y0 = 0;
z0 = 0;
vx0 = 0;
vy0 = 0.8;
vz0 = 0;

r0 = sqrt((x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2);
Er = -GM*m/r0 + 0.5*m*(vx0^2 + vy0^2 + vz0^2);
% a = 1; % simulate the circular orbit first
RE = -GM*m/2/Er;
vRE = sqrt(2*GM/RE - GM/RE);


% Set the angle of rotation to be anti-clockwise viewing from positive z to negative z direction.

%% use Newton's law to calculate locations and speed
%initial condition
dataN = zeros(length(t), 9);  % t, x, y, vx, vy, r, Er
dataN(:,1) = t; 
dataN(1, 2) = x0;
dataN(1, 3) = y0;
dataN(1, 4) = z0;
dataN(1, 5) = vx0;
dataN(1, 6) = vy0;
dataN(1, 7) = vz0;
dataN(1, 8) = RE;
dataN(1, 9) = Er;

% calculate the state stepwisely
for i = 2:length(t)
    %read previous data
    x1 = dataN(i-1, 2); y1 = dataN(i-1, 3); z1 = dataN(i-1, 4);
    vx1 = dataN(i-1, 5); vy1 = dataN(i-1, 6); vz1 = dataN(i-1, 7);
    r1 = dataN(i-1, 8); Er1 = dataN(i-1, 9);
    [Fx, Fy, Fz] = force(x1, y1, z1, xc, yc, zc, GM, m);
    ax = Fx/m; ay = Fy/m; az = Fz/m;
    vx2 = vx1 + ax*dt; vy2 = vy1 + ay*dt; vz2 = vz1 + az*dt;
    x2 = x1 + vx2*dt; y2 = y1 + vy2*dt; z2 = z1 + vz2*dt;
    
    r = sqrt((x2-xc)^2 + (y2-yc)^2 + (z2-zc)^2);
    E = -GM*m/r + 0.5*m*(vx2^2 + vy2^2 + vz2^2);
    
    dataN(i, 2) = x2;
    dataN(i, 3) = y2;
    dataN(i, 4) = z2;
    dataN(i, 5) = vx2;
    dataN(i, 6) = vy2;
    dataN(i, 7) = vz2;
    dataN(i, 8) = r;
    dataN(i, 9) = E;    
end

% figure; plot3(dataN(:,2), dataN(:,3), dataN(:,4));

figure;
plot(dataN(:,2), dataN(:,3), 'r-'); hold on;
plot(xc, yc, 'k+', 'markersize', 10, 'linewidth', 2);


%% use vis-viva equation and construct a complex number for the state
% initialize parameters
GM = 1;
m = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be origin and stationary

%initialize data
dt = 0.001;
t = 0:dt:6;
t = t';
x0 = 1;
y0 = 0;
z0 = 0;
vx0 = 0;   %only works for starting from apogee and in x-y plane anticlockwise now. Need revision to work for random location and angle.
vy0 = 0.8;  
vz0 = 0;
r0 = sqrt( (x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2 );
Er0 = -GM*m/r0 + 0.5*m*(vx0^2 + vy0^2 + vz0^2);
RE = -GM*m/2/Er0;

vcr = cross([x0-xc, y0-yc, z0-zc], [vx0, vy0, z0-zc]);
vcr = sqrt(dot(vcr, vcr));
a0 = vcr^2/GM/RE;

dataV = zeros(length(t), 6);  % t, x, y, z, wf, r, Er
dataV(:,1) = t;
dataV(1, 2) = x0;
dataV(1, 3) = y0;
dataV(1, 4) = findwf(x0, y0, xc, yc, Er0, a0, GM, m);
dataV(1, 5) = r0;
dataV(1, 6) = Er0;


%calculate state stepwisely
for i = 2:length(t)
    x1 = dataV(i-1, 2); y1 = dataV(i-1, 3);
    phi1 = dataV(i-1, 4);
    r1 = dataV(i-1, 5); Er1 = dataV(i-1, 6);
    
    theta1 = atan2d(y1-yc, x1-xc)*pi()/180;
    thetad = real(phi1)/m*dt/r1;
    theta2 = theta1 + thetad;
    if y1 >= 0
        approach = true;
    else
        approach = false;
    end
    if approach
        r2 = r1*cos(thetad) - imag(phi1)/m*dt; % - OR + ???
    else
        r2 = r1/cos(thetad) + imag(phi1)/m*dt; % - OR + ???
    end
    x2 = r2*cos(theta2) + xc;
    y2 = r2*sin(theta2) + yc;

    dataV(i, 2) = x2;
    dataV(i, 3) = y2;
    dataV(i, 4) = findwf(x2, y2, xc, yc, Er0, a0, GM, m);
    dataV(i, 5) = sqrt((x2-xc)^2 + (y2-yc)^2);
    dataV(i, 6) = Er0;
end

figure; 
plot(dataV(:,2), dataV(:,3), 'b-'); hold on;
plot(xc, yc, 'k+', 'markersize', 10, 'linewidth', 2);


%% make the movie

figure; 
plot(dataV(:,2), dataV(:,3), 'b-'); hold on;
plot(xc, yc, 'k+', 'markersize', 10, 'linewidth', 2);

for i = 1:length(t)
    plot(dataN(i,2), dataN(i,3), 'rd');
    plot(dataV(i,2), dataV(i,3), 'bo');
    jcPlotStyle;
    pause(0.001)
end
hold off



% Video = VideoWriter('one-body.avi');
% open(Video);
% 
% 
% h = figure; 
% plot(dataV(:,2), dataV(:,3), 'b-'); hold on;
% plot(xc, yc, 'k+', 'markersize', 10, 'linewidth', 2);
% 
% for i = 1:10:length(t)
%     plot(dataN(i,2), dataN(i,3), 'rd');
%     plot(dataV(i,2), dataV(i,3), 'bo');
%     jcPlotStyle;
%     pause(0.01)
%     Frame = getframe(h);
%     writeVideo(Video, Frame);
% end
% hold off



%% functions
% find the wavefunction:
function phi = findwf(x1, y1, xc, yc, Er0, a0, GM, m)
    r = sqrt((x1-xc)^2 + (y1-yc)^2);
    RE = -GM*m/2/Er0;
    vp = (sqrt(a0*GM*RE)/r);
    vr = (sqrt(-GM/RE + 2*GM/r - a0*GM*RE/r^2));
    phi = m*vp + 1i*m*vr;
end


% function to calculate the gravitational force
function [Fx, Fy, Fz] = force(x, y, z, xc, yc, zc, GM, m)
   rx = x - xc;
   ry = y - yc;
   rz = z - zc;
   r = sqrt(rx^2 + ry^2 + rz^2);
   F = GM*m/r^2;
   Fx = -rx/r*F;
   Fy = -ry/r*F;    
   Fz = -rz/r*F;
end



