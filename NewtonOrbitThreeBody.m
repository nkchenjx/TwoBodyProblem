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
% close all;

% initialize parameters
G = 1;
m1 = 1;
m2 = 1;
m3 = 1;
xc = 0; yc = 0; zc = 0; %center of mass set to be stationary as the inertia reference frame, thus two bodies have correlated freedom.

dt = 0.001;
t = 0:dt:10;
t = t';

% %initialize data -- stable solution of number 8
% x10 = -0.97000436;  % initialize state of m1
% y10 = 0.24308753;
% z10 = 0;
% vx10 = 0.4662036850;
% vy10 = 0.4323657300;
% vz10 = 0;
% x20 = 0;  % initialize state of  m2
% y20 = 0;
% z20 = 0;
% vx20 = -0.93240737;
% vy20 = -0.86473146;
% vz20 = 0;

%initialize data -- stable solution of 
x10 = -0.97000436;  % initialize state of m1
y10 = 0.24308753;
z10 = 0;
vx10 = 0.4662036850;
vy10 = 0.4323657300;
vz10 = 0;
x20 = 0;  % initialize state of  m2
y20 = 0;
z20 = 0;
vx20 = -0.93240737;
vy20 = -0.86473146;
vz20 = 0.2;


x30 = (xc-x10)*m1/m3 + (xc-x20)*m2/m3 + xc; % state of m3 is clculatied to keep the center of mass stationary
y30 = (yc-y10)*m1/m3 + (yc-y20)*m2/m3 + yc;
z30 = (zc-z10)*m1/m3 + (zc-z20)*m2/m3 + zc;
vx30 = -vx10*(m1/m3) - vx20*(m2/m3);
vy30 = -vy10*(m1/m3) - vy20*(m2/m3);
vz30 = -vz10*(m1/m3) - vz20*(m2/m3);


r120 = sqrt((x10-x20)^2 + (y10-y20)^2 + (z10-z20)^2);
r130 = sqrt((x10-x30)^2 + (y10-y30)^2 + (z10-z30)^2);
r230 = sqrt((x20-x30)^2 + (y20-y30)^2 + (z20-z30)^2);

figure; hold on; 
plot(xc, yc, '+', 'markersize', 10, 'linewidth', 2); 
plot(x10, y10, 'o', 'markersize', 10, 'linewidth', 2); 
plot(x20, y20, 'o', 'markersize', 10, 'linewidth', 2); 
plot(x30, y30, 'o', 'markersize', 10, 'linewidth', 2);
hold off;

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

data3 = zeros(length(t), 7);  % t, x, y, z, vx, vy, vz
data3(:,1) = t; 
data3(1, 2) = x30;
data3(1, 3) = y30;
data3(1, 4) = z30;
data3(1, 5) = vx30;
data3(1, 6) = vy30;
data3(1, 7) = vz30;



dataC = zeros(length(t), 4); %t, Center of mass (CoM) x, y, z, and and miu23, miu13, miu12
dataC(:,1) = t; 
dataC(1, 2) = xc + (x10*m1 + x20*m2 + x30*m3)/(m1 + m2 + m3); %CoM
dataC(1, 3) = yc + (y10*m1 + y20*m2 + y30*m3)/(m1 + m2 + m3);
dataC(1, 4) = zc + (z10*m1 + z20*m2 + z30*m3)/(m1 + m2 + m3);
i = 1;
dataC(i,5:13) = findmiu(m1, m2, m3, xc, yc, zc, data1(i,2:4), data2(i,2:4), data3(i,2:4));



% calculate the state stepwisely
for i = 2:length(t)
    %read previous data
    x11 = data1(i-1, 2); y11 = data1(i-1, 3); z11 = data1(i-1, 4);
    vx11 = data1(i-1, 5); vy11 = data1(i-1, 6); vz11 = data1(i-1, 7);
    x21 = data2(i-1, 2); y21 = data2(i-1, 3); z21 = data2(i-1, 4);
    vx21 = data2(i-1, 5); vy21 = data2(i-1, 6); vz21 = data2(i-1, 7);
    x31 = data3(i-1, 2); y31 = data3(i-1, 3); z31 = data3(i-1, 4);
    vx31 = data3(i-1, 5); vy31 = data3(i-1, 6); vz31 = data3(i-1, 7);
    
    r12 = sqrt((x11-x21)^2 + (y11-y21)^2 + (z11-z21)^2);
    r13 = sqrt((x11-x31)^2 + (y11-y31)^2 + (z11-z31)^2);
    r23 = sqrt((x21-x31)^2 + (y21-y31)^2 + (z21-z31)^2);
    if r12
        Fx12 = G*m1*m2*(x21-x11)/r12^3; Fx21 = G*m1*m2*(x11-x21)/r12^3;
        Fy12 = G*m1*m2*(y21-y11)/r12^3; Fy21 = G*m1*m2*(y11-y21)/r12^3;
        Fz12 = G*m1*m2*(z21-z11)/r12^3; Fz21 = G*m1*m2*(z11-z21)/r12^3;
    else %r12 = 0, singularity
        Fx12 = 0; Fx21 = 0;
        Fy12 = 0; Fy21 = 0;
        Fz12 = 0; Fz21 = 0;
    end
    if r13
        Fx13 = G*m1*m3*(x31-x11)/r13^3; Fx31 = G*m1*m3*(x11-x31)/r13^3;
        Fy13 = G*m1*m3*(y31-y11)/r13^3; Fy31 = G*m1*m3*(y11-y31)/r13^3;
        Fz13 = G*m1*m3*(z31-z11)/r13^3; Fz31 = G*m1*m3*(z11-z31)/r13^3;
    else %r13 = 0, singularity
        Fx13 = 0; Fx31 = 0;
        Fy13 = 0; Fy31 = 0;
        Fz13 = 0; Fz31 = 0;
    end
    if r23
        Fx23 = G*m3*m2*(x31-x21)/r23^3; Fx32 = G*m3*m2*(x21-x31)/r23^3;
        Fy23 = G*m3*m2*(y31-y21)/r23^3; Fy32 = G*m3*m2*(y21-y31)/r23^3;
        Fz23 = G*m3*m2*(z31-z21)/r23^3; Fz32 = G*m3*m2*(z21-z31)/r23^3;
    else %r23 = 0, singularity
        Fx23 = 0; Fx32 = 0;
        Fy23 = 0; Fy32 = 0;
        Fz23 = 0; Fz32 = 0;
    end  
    
    
    ax1 = (Fx12 + Fx13)/m1; ay1 = (Fy12 + Fy13)/m1; az1 = (Fz12 + Fz13)/m1;
    ax2 = (Fx21 + Fx23)/m2; ay2 = (Fy21 + Fy23)/m2; az2 = (Fz21 + Fz23)/m2;
    ax3 = (Fx31 + Fx32)/m3; ay3 = (Fy31 + Fy32)/m3; az3 = (Fz31 + Fz32)/m3;
    
    vx12 = vx11 + ax1*dt; vy12 = vy11 + ay1*dt; vz12 = vz11 + az1*dt;
    vx22 = vx21 + ax2*dt; vy22 = vy21 + ay2*dt; vz22 = vz21 + az2*dt;
    vx32 = vx31 + ax3*dt; vy32 = vy31 + ay3*dt; vz32 = vz31 + az3*dt;
    
    x12 = x11 + vx12*dt; y12 = y11 + vy12*dt; z12 = z11 + vz12*dt;
    x22 = x21 + vx22*dt; y22 = y21 + vy22*dt; z22 = z21 + vz22*dt;
    x32 = x31 + vx32*dt; y32 = y31 + vy32*dt; z32 = z31 + vz32*dt;    
    
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
    
    data3(i, 2) = x32; % m3 state
    data3(i, 3) = y32;
    data3(i, 4) = z32;
    data3(i, 5) = vx32;
    data3(i, 6) = vy32;
    data3(i, 7) = vz32;
    

    dataC(i, 2) = xc + (x12*m1 + x22*m2 + x32*m3)/(m1 + m2 + m3); %CoM %CoM should not change but check for error.
    dataC(i, 3) = yc + (y12*m1 + y22*m2 + y32*m3)/(m1 + m2 + m3);
    dataC(i, 4) = zc + (z12*m1 + z22*m2 + z32*m3)/(m1 + m2 + m3);
    dataC(i,5:13) = findmiu(m1, m2, m3, xc, yc, zc, data1(i,2:4), data2(i,2:4), data3(i,2:4));
end

figure; hold on
plot3(data1(:,2), data1(:,3), data1(:,4));
plot3(data2(:,2), data2(:,3), data2(:,4));
plot3(data3(:,2), data3(:,3), data3(:,4));
plot3(dataC(:,2), dataC(:,3), dataC(:,4), 'k+', 'markersize', 10, 'linewidth', 2); %CoM
hold off
% plot3(dataC(:,5), dataC(:,6), dataC(:,7),'kx'); %CoE
figure; plot(dataC(:,1), dataC(:,5:13));

pause(1);

%% calculate angular momentum
data = cat(3, data1, data2, data3);
L = [];
for i = 1:size(data, 1)
    L(i,:) = findAngularMomentum([m1, m2, m3], xc, yc, zc, data(i, :, :)); 
end

figure; plot(L(:,1), L(:,2:end));
pause(1);



%% make the movie

figure; hold on
plot3(data1(:,2), data1(:,3), data1(:,4), 'r');
plot3(data2(:,2), data2(:,3), data2(:,4),'b');
plot3(data3(:,2), data3(:,3), data3(:,4),'g');
plot3(dataC(:,2), dataC(:,3), dataC(:,4), 'k+', 'markersize', 10, 'linewidth', 2); %CoM

for i = 1:10:length(t)
    plot3(data1(i,2), data1(i,3), data1(i,4), 'ro');
    plot3(data2(i,2), data2(i,3), data2(i,4), 'bo');
    plot3(data3(i,2), data3(i,3), data3(i,4), 'go');
    plot3(dataC(i,2), dataC(i,3), dataC(i,4), 'k+', 'markersize', 10, 'linewidth', 2); %CoM
%     plot(dataC(i,5), dataC(i,6), 'kx'); %CoE 
    pause(0.01);
end
hold off



%% functions
function mius = findmiu(m1, m2, m3, xc, yc, zc, data1, data2, data3) % miu 23, miu13, miu12
    x1 = data1(1); y1 = data1(2); z1 = data1(3);
    x2 = data2(1); y2 = data2(2); z2 = data2(3);
    x3 = data3(1); y3 = data3(2); z3 = data3(3);
    r12 = sqrt(dot(data1 - data2, data1 - data2));
    r21 = r12;
    r13 = sqrt(dot(data1 - data3, data1 - data3));
    r31 = r13;
    r23 = sqrt(dot(data2 - data3, data2 - data3));
    r32 = r23;
    r1 = sqrt(dot(data1, data1));
    r2 = sqrt(dot(data2, data2));
    r3 = sqrt(dot(data3, data3));
    
    mius(1, 1) = (m2*r1^3*(x2-x1)/r12^3 + m3*(x3-x1)*r1^3/r13^3)/(xc-x1); %miu23x
    mius(1, 2) = (m1*r2^3*(x1-x2)/r21^3 + m3*(x3-x2)*r2^3/r23^3)/(xc-x2); %miu13x
    mius(1, 3) = (m1*r3^3*(x1-x3)/r13^3 + m2*(x2-x3)*r3^3/r32^3)/(xc-x3); %miu12x
    
    mius(1, 4) = (m2*r1^3*(y2-y1)/r12^3 + m3*(y3-y1)*r1^3/r13^3)/(yc-y1); %miu23y
    mius(1, 5) = (m1*r2^3*(y1-y2)/r21^3 + m3*(y3-y2)*r2^3/r23^3)/(yc-y2); %miu13y
    mius(1, 6) = (m1*r3^3*(y1-y3)/r13^3 + m2*(y2-y3)*r3^3/r32^3)/(yc-y3); %miu12y

    mius(1, 7) = (m2*r1^3*(z2-z1)/r12^3 + m3*(z3-z1)*r1^3/r13^3)/(zc-z1); %miu23z
    mius(1, 8) = (m1*r2^3*(z1-z2)/r21^3 + m3*(z3-z2)*r2^3/r23^3)/(zc-z2); %miu13z
    mius(1, 9) = (m1*r3^3*(z1-z3)/r13^3 + m2*(z2-z3)*r3^3/r32^3)/(zc-z3); %miu12z

 
    
end

function L = findAngularMomentum(m, xc, yc, zc, data) 
    %data: t, x, y, z, vx, vy, vz
    % m: [m1, m2, ...]
    numP = size(data, 3);
    Lt = [0,0,0];
    for i = 1:numP
        x = data(1, 2, i) - xc;
        y = data(1, 3, i) - xc;
        z = data(1, 4, i) - xc;
        vx = data(1, 5, i);
        vy = data(1, 6, i);
        vz = data(1, 7, i);
        
        L(1, (i-1)*3+1) = m(i)*(y*vz - z*vy); % Lx
        L(1, (i-1)*3+2) = m(i)*(z*vx - x*vz); % Ly
        L(1, (i-1)*3+3) = m(i)*(x*vy - y*vx); % Lz
        Lt(1) = Lt(1) + L(1, (i-1)*3+1);
        Lt(2) = Lt(2) + L(1, (i-1)*3+2);
        Lt(3) = Lt(3) + L(1, (i-1)*3+3);
    end
    
    L = [data(1), L, Lt];
        
end
