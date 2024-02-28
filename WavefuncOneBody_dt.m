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

clear;
close all;


%% use vis-viva equation and construct a complex number for the state
% initialize parameters
GM = 1;
m = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be origin and stationary

%initialize data
dt = 0.001;
t = 0:dt:20;
t = t';
x0 = 1;
y0 = 0;
z0 = 0;
vx0 = -0.2;   %only works for starting from apogee and in x-y plane anticlockwise now. Need revision to work for random location and angle.
vy0 = 0.6;  
vz0 = 0;
r0 = sqrt( (x0-xc)^2 + (y0-yc)^2 + (z0-zc)^2 );
Er0 = -GM*m/r0 + 0.5*m*(vx0^2 + vy0^2 + vz0^2);
RE = -GM*m/2/Er0;

vcr = cross([x0-xc, y0-yc, z0-zc], [vx0, vy0, vz0]);
vcr = sqrt(dot(vcr, vcr));
a0 = vcr^2/GM/RE;

dataV = zeros(length(t), 6);  % t, x, y, wf, r, Er, theta
dataV(:,1) = t;
dataV(1, 2) = x0;
dataV(1, 3) = y0;
dataV(1, 4) = r0;
dataV(1, 5) = Er0;
dataV(1, 6) = 0; %theta
phi = zeros(length(t), 1);


if vx0 <= 0
    sign = -1; %aproaching center
else
    sign = 1;
end
phi(1, 1) = findwf(x0, y0, xc, yc, Er0, a0, GM, m, sign); % using complex number slows down calculation so maybe split it to a different array or use mvr mvi instead.

[rapogee, rperigee] = findsemiaxis(x0, y0, z0, xc, yc, zc, vx0, vy0, vz0, GM, m);


%calculate state stepwisely
for i = 2:length(t)
    x1 = dataV(i-1, 2); y1 = dataV(i-1, 3);
    phi1 = phi(i-1);
    r1 = dataV(i-1, 4); Er1 = dataV(i-1, 5);
    
    theta1 = atan2d(y1-yc, x1-xc)*pi()/180;
    thetad = real(phi1)*dt/m/r1;
    theta2 = theta1 + thetad;
    
    if r1 <= rperigee 
        sign = 1;
    end
    if r1 >= rapogee 
        sign = -1;
    end

    r2 = r1 + sign*abs(imag(phi1))/m*dt; % - OR + ???

    x2 = r2*cos(theta2) + xc;
    y2 = r2*sin(theta2) + yc;

    dataV(i, 2) = x2;
    dataV(i, 3) = y2;
    dataV(i, 4) = r2;
    dataV(i, 5) = Er0;
    dataV(i, 6) = theta2;
    phi(i) = findwf(x2, y2, xc, yc, Er0, a0, GM, m, sign);
end

figure; 
plot(dataV(:,2), dataV(:,3)); hold on;
plot(xc, yc, 'k+', 'markersize', 10, 'linewidth', 2); 

figure; plot(real(dataV(:, 4)), real(phi(:))); hold on; 
plot(real(dataV(:, 4)), imag(phi(:))); xlabel('r');

figure; plot(real(dataV(:, 1)), real(phi(:))); hold on; 
plot(real(dataV(:, 1)), imag(phi(:))); xlabel('t');

figure; plot(real(dataV(:, 6)), real(phi(:))); hold on; 
plot(real(dataV(:, 6)), imag(phi(:))); xlabel('theta');

%% functions
% find the wavefunction:
function phi = findwf(x1, y1, xc, yc, Er0, a0, GM, m, sign)
    r = sqrt((x1-xc)^2 + (y1-yc)^2);
    RE = -GM*m/2/Er0;
    vc = abs(sqrt(a0*GM*RE)/r); %some case near 0 has a leak and cause v to be imaginary, so add abs.
    vr = abs(sqrt(-GM/RE + 2*GM/r - a0*GM*RE/r^2));

    phi = m*vc + sign*1i*m*vr;

end
 
function [ra, rp] = findsemiaxis(x, y, z, xc, yc, zc, vx, vy, vz, GM, m)
    r = sqrt( (x-xc)^2 + (y-yc)^2 + (z-zc)^2 );
    Er = -GM*m/r + 0.5*m*(vx^2 + vy^2 + vz^2);
    RE = -GM*m/2/Er;
    
    vcr = cross([x-xc, y-yc, z-zc], [vx, vy, vz]);
    vcr = sqrt(dot(vcr, vcr));
    a0 = vcr^2/GM/RE;
    
    a = -a0*GM*RE;
    b = 2*GM;
    c = -GM/RE;
    
    ra = 2*a / (-b + sqrt(b^2-4*a*c));
    rp = 2*a / (-b - sqrt(b^2-4*a*c));
end