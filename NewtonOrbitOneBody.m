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

% clear;

% initialize parameters
GM = 1;
m = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be origin and stationary

%initialize data
dt = 0.001;
t = 0:dt:2.75;
t = t';
x0 = 1;
y0 = 0;
z0 = 0;
vx0 = 0;
vy0 = 0.3;
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

figure; plot(dataN(:,2), dataN(:,3));

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

