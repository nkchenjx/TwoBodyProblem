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

% Caution only x-y plane motion is simulated, z has to be zero. also only apogee at
% x-axis is simulated. Random orientation and starting point have not been
% coded.

% clear;
% close all;


%% use vis-viva equation and construct a complex number for the state
% initialize parameters
GM = 1;
m = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be origin and stationary

%initialize data
dt = 0.0001;
t = 0:dt:5;
t = t';
x0 = 1;
y0 = 0;
z0 = 0;
vx0 = 0;   %only works for starting from apogee and in x-y plane anticlockwise now. Need revision to work for random location and angle.
vy0 = 0.6;  
vz0 = 0;
r0 = sqrt( (x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2 );
Er0 = -GM*m/r0 + 0.5*m*(vx0^2 + vy0^2 + vz0^2);
RE = -GM*m/2/Er0;

vcr = cross([x0-xc, y0-yc, z0-zc], [vx0, vy0, vz0]);
vcr = sqrt(dot(vcr, vcr));
a0 = vcr^2/GM/RE;

dataV = zeros(length(t), 7);  % t, x, y, z, wf, r, Er, theta
dataV(:,1) = t;
dataV(1, 2) = x0;
dataV(1, 3) = y0;
dataV(1, 4) = findwf(x0, y0, xc, yc, Er0, a0, GM, m);
dataV(1, 5) = r0;
dataV(1, 6) = Er0;
dataV(1, 7) = 0;


%calculate state stepwisely
for i = 2:length(t)
    x1 = dataV(i-1, 2); y1 = dataV(i-1, 3);
    phi1 = dataV(i-1, 4);
    r1 = dataV(i-1, 5); Er1 = dataV(i-1, 6);
    
    theta1 = atan2d(y1-yc, x1-xc)*pi()/180;
    thetad = real(phi1)/m*dt/r1;
    theta2 = theta1 + thetad;

    r2 = r1/cos(thetad) + imag(phi1)/m*dt; % - OR + ???
 
    x2 = r2*cos(theta2) + xc;
    y2 = r2*sin(theta2) + yc;

    dataV(i, 2) = x2;
    dataV(i, 3) = y2;
    dataV(i, 4) = findwf(x2, y2, xc, yc, Er0, a0, GM, m);
    dataV(i, 5) = sqrt((x2-xc)^2 + (y2-yc)^2);
    dataV(i, 6) = Er0;
    dataV(i, 7) = theta2;
end

figure; 
plot(dataV(:,2), dataV(:,3));

figure; plot(real(dataV(:, 5)), real(dataV(:,4))); hold on; 
plot(real(dataV(:, 5)), imag(dataV(:,4))); xlabel('r');

figure; plot(real(dataV(:, 1)), real(dataV(:,4))); hold on; 
plot(real(dataV(:, 1)), imag(dataV(:,4))); xlabel('t');

figure; plot(real(dataV(:, 7)), real(dataV(:,4))); hold on; 
plot(real(dataV(:, 7)), imag(dataV(:,4))); xlabel('theta');

%% functions
% find the wavefunction:
function phi = findwf(x1, y1, xc, yc, Er0, a0, GM, m)
    r = sqrt((x1-xc)^2 + (y1-yc)^2);
    RE = -GM*m/2/Er0;
    vp = abs(sqrt(a0*GM*RE)/r); %some case near 0 has a leak and cause v to be imaginary, so add abs.
    vr = abs(sqrt(-GM/RE + 2*GM/r - a0*GM*RE/r^2));
    if y1 >= 0  % only works if x axis is the long semi-axis and apogee is on the right. Need transformation for other cases.
        approach = true; % moving away from apogee
        phi = m*vp - 1i*m*vr;
    else
        approach = false;
        phi = m*vp + 1i*m*vr;
    end
 end