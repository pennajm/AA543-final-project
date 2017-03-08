%AA543 HW3 grid generation
clear all;
close all;

afbody = importdata('body.dat');
Imax = 129;
Jmax = 65;
Is = 1:1:Imax;
Js = 1:1:Jmax;

%pre-allocate arrays
x = zeros(Imax, Jmax);
y = zeros(Imax, Jmax);
Ax = 0.0; %zeros(Imax, Jmax);
Ay = 0.0;%zeros(Imax, Jmax);
Bx = 0.0;%zeros(Imax, Jmax);
By = 0.0;%zeros(Imax, Jmax);
Tx = 0.0;%zeros(Imax, Jmax);
Ty = 0.0;%zeros(Imax, Jmax);
L1x = 0.0;%zeros(1, Imax);
L2x = 0.0;%zeros(1, Imax);
L1y = 0.0;%zeros(1, Jmax);
L2y = 0.0;%zeros(1, Jmax);
x(:,1) = afbody(:,2) -0.5;
y(:,1) = afbody(:,3);
%set up outer boundary
clen= 1.0;
% circle
r=10*clen; % radius
theta0=0:2*pi/(Imax-1):2*pi; % the angle
theta=linspace(2*pi, 0, 129);
%theta = flip(theta0);
x(:,Jmax) = r*cos(theta');
y(:, Jmax) = r*sin(theta');
%now the chord cuts
x(1,:) = linspace(0.5*clen, (10*clen), Jmax);
y(1, :) = zeros(1,Jmax);
x(Imax, :) = linspace(0.5*clen, (10*clen), Jmax);
y(Imax, :) = zeros(1,Jmax);
x0 = x;
y0 = y;
%plot boundaries
%figure(1);
%plot(x0,y0);
%set up the initial conditions with transfinite interpolation
for jj=2:(Jmax-1)
    for ii = 2:(Imax-1)
        L1x = (ii -Imax)/(1 - Imax);
        L2x = (ii -1)/(Imax -1);
        L1y =  (jj - Jmax)/(1 - Jmax);
        L2y = (jj -1)/(Jmax -1);
        Ax = L1x*x(1, jj) + L2x*x(Imax, jj);
        Ay = L1x*y(1, jj) + L2x*y(Imax, jj);
        Bx = L1y*x(ii, 1) + L2y*x(ii, Jmax);
        By = L1y*y(ii, 1) + L2y*y(ii, Jmax);
        Tx = L1x*L1y*x(1, 1) + L2x*L2y*x(Imax, Jmax) ...
            + L1x*L2y*x(1, Jmax) + L2x*L1y*x(Imax,1);
        Ty = L1x*L1y*y(1, 1) + L2x*L2y*y(Imax, Jmax) ...
            + L1x*L2y*y(1, Jmax) + L2x*L1y*y(Imax,1);
        x(ii,jj) = Ax + Bx -Tx;
        y(ii,jj) = Ay + By -Ty;
    end
end
%figure(1);
%plot(x,y);
%initial conditions
xi = x;
yi = y;
%set omega for SOR
w = 1.25;
%coefficients for the P and Q
%first set
PQ1 = [0,0]; %trivial case where P =Q =0
%second set
%first set coeffs
a = 15.0;%linspace(-10, 10, Imax);
b = 1.0;%linspace(-10,10, Jmax);
c = 1.0;%linspace(1,10, Imax);
d = 1.0;%linspace(1,10, Jmax);
etam = 1;%ones(1, Imax);
etal = 1;
xil = 1;
xim = 1;
PQ2 = zeros(2, Imax, Jmax);
for jj=1:Jmax
    for ii=1:Imax
        PQ2(1,ii,jj) = -sum(a.*sign(ii - xil).*exp(-c.*abs(ii -xil))) ...
        -sum(b.*sign(ii - xil).*exp(-d.*sqrt((ii -xim).^2 + (jj -etam).^2)));

        PQ2(2,ii,jj) = -sum(a.*sign(jj - etal).*exp(-c.*abs(jj -etal))) ...
        -sum(b.*sign(jj - etam).*exp(-d.*sqrt((ii -xim).^2 + (jj -etam).^2)));
    end
end
%now employ successive overrelaxation
n =0;
R = 0;
%Main loop
ResAr = [];
Res = 1.0;
while Res > 1e-6
    n = n +1;
    for jj = 2:(Jmax-1)
        for ii = 1:Imax
            xold = x(ii, jj);
            yold = y(ii, jj);
            %P = PQ1(1);
            %Q = PQ1(2);
            P = PQ2(1, ii, jj);
            Q = PQ2(2, ii, jj);
            if ii ==1
                ap = 0.25*((x(ii, jj+1) - x(ii, jj-1))^2 + (y(ii, jj+1) - y(ii, jj-1))^2);
                bp = 0.25*((x(ii+1, jj) - x(Imax-1, jj))*(x(ii, jj+1) - x(ii, jj-1)) ...
                    + (y(ii+1, jj) - y(Imax-1, jj))*(y(ii, jj+1) - y(ii, jj-1)));
                gp = 0.25*((x(ii+1, jj) - x(Imax-1, jj))^2 + (y(ii+1, jj) - y(Imax-1, jj))^2);
                dp = 0.0625*((x(ii+1, jj) - x(Imax-1, jj))*(y(ii, jj+1) - y(ii, jj-1)) ...
                    - (x(ii, jj+1) - x(ii, jj-1))*(y(ii+1, jj) - y(Imax-1, jj)))^2;
                fcoeff = 1.0/(2.0*ap + 2.0*gp);
                xbarij = fcoeff*(ap*(x(Imax-1, jj) +x(ii+1,jj)) - 0.5*bp*(x(ii+1, jj+1) ...
                    -x(Imax-1, jj+1) - x(ii+1, jj-1) + x(Imax-1, jj-1)) + gp*(x(ii, jj-1) + x(ii, jj+1)) ...
                    +0.5*dp*(P*(x(ii+1, jj) - x(Imax-1,jj)) + Q*(x(ii, jj+1) - x(ii, jj-1))));
                ybarij = fcoeff*(ap*(y(Imax-1, jj) +y(ii+1,jj)) - 0.5*bp*(y(ii+1, jj+1) ...
                    -y(Imax-1, jj+1) - y(ii+1, jj-1) + y(Imax-1, jj-1)) + gp*(y(ii, jj-1) + y(ii, jj+1)) ...
                    +0.5*dp*(P*(y(ii+1, jj) - y(Imax-1,jj)) + Q*(y(ii, jj+1) - y(ii, jj-1))));
            elseif ii ==Imax
                ap = 0.25*((x(ii, jj+1) - x(ii, jj-1))^2 + (y(ii, jj+1) - y(ii, jj-1))^2);
                bp = 0.25*((x(2, jj) - x(ii-1, jj))*(x(ii, jj+1) - x(ii, jj-1)) ...
                    + (y(2, jj) - y(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)));
                gp = 0.25*((x(2, jj) - x(ii-1, jj))^2 + (y(2, jj) - y(ii-1, jj))^2);
                dp = 0.0625*((x(2, jj) - x(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)) ...
                    - (x(ii, jj+1) - x(ii, jj-1))*(y(2, jj) - y(ii-1, jj)))^2;
                fcoeff = 1.0/(2.0*ap + 2.0*gp);
                xbarij = fcoeff*(ap*(x(ii-1, jj) +x(2,jj)) - 0.5*bp*(x(2, jj+1) ...
                    -x(ii-1, jj+1) - x(2, jj-1) + x(ii-1, jj-1)) + gp*(x(ii, jj-1) + x(ii, jj+1)) ...
                    +0.5*dp*(P*(x(2, jj) - x(ii-1,jj)) + Q*(x(ii, jj+1) - x(ii, jj-1))));
                ybarij = fcoeff*(ap*(y(ii-1, jj) +y(2,jj)) - 0.5*bp*(y(2, jj+1) ...
                    -y(ii-1, jj+1) - y(2, jj-1) + y(ii-1, jj-1)) + gp*(y(ii, jj-1) + y(ii, jj+1)) ...
                    +0.5*dp*(P*(y(2, jj) - y(ii-1,jj)) + Q*(y(ii, jj+1) - y(ii, jj-1))));
            else
                %coefficients
                ap = 0.25*((x(ii, jj+1) - x(ii, jj-1))^2 + (y(ii, jj+1) - y(ii, jj-1))^2);
                bp = 0.25*((x(ii+1, jj) - x(ii-1, jj))*(x(ii, jj+1) - x(ii, jj-1)) ...
                    + (y(ii+1, jj) - y(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)));
                gp = 0.25*((x(ii+1, jj) - x(ii-1, jj))^2 + (y(ii+1, jj) - y(ii-1, jj))^2);
                dp = 0.0625*((x(ii+1, jj) - x(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)) ...
                    - (x(ii, jj+1) - x(ii, jj-1))*(y(ii+1, jj) - y(ii-1, jj)))^2;
                fcoeff = 1.0/(2.0*ap + 2.0*gp);
                xbarij = fcoeff*(ap*(x(ii-1, jj) +x(ii+1,jj)) - 0.5*bp*(x(ii+1, jj+1) ...
                    -x(ii-1, jj+1) - x(ii+1, jj-1) + x(ii-1, jj-1)) + gp*(x(ii, jj-1) + x(ii, jj+1)) ...
                    +0.5*dp*(P*(x(ii+1, jj) - x(ii-1,jj)) + Q*(x(ii, jj+1) - x(ii, jj-1))));
                ybarij = fcoeff*(ap*(y(ii-1, jj) +y(ii+1,jj)) - 0.5*bp*(y(ii+1, jj+1) ...
                    -y(ii-1, jj+1) - y(ii+1, jj-1) + y(ii-1, jj-1)) + gp*(y(ii, jj-1) + y(ii, jj+1)) ...
                    +0.5*dp*(P*(y(ii+1, jj) - y(ii-1,jj)) + Q*(y(ii, jj+1) - y(ii, jj-1))));
            end
            x(ii,jj) = w*xbarij + (1-w)*xold;
            y(ii,jj) = w*ybarij + (1-w)*yold;
        end
    end
    R = 0.0;
    for jj=2:(Jmax-1)
        for ii=2:(Imax-1)
            ap = 0.25*((x(ii, jj+1) - x(ii, jj-1))^2 + (y(ii, jj+1) - y(ii, jj-1))^2);
            bp = 0.25*((x(ii+1, jj) - x(ii-1, jj))*(x(ii, jj+1) - x(ii, jj-1)) ...
                + (y(ii+1, jj) - y(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)));
            gp = 0.25*((x(ii+1, jj) - x(ii-1, jj))^2 + (y(ii+1, jj) - y(ii-1, jj))^2);
            dp = 0.0625*((x(ii+1, jj) - x(ii-1, jj))*(y(ii, jj+1) - y(ii, jj-1)) ...
                + (x(ii, jj+1) - x(ii, jj-1))*(y(ii+1, jj) - y(ii-1, jj)))^2;
            P = PQ1(1);
            Q = PQ1(2);
            Rx= ap*(x(ii-1, jj) -2*x(ii, jj) + x(ii+1, jj)) ...
                -0.5*bp*(x(ii+1, jj+1) -x(ii-1, jj+1) - x(ii+1, jj-1) + x(ii-1, jj-1)) ...
                +gp*(x(ii, jj-1) -2*x(ii, jj) + x(ii, jj+1)) + 0.5*dp*(P*(x(ii+1, jj) - x(ii-1,jj)) ...
                + Q*(x(ii, jj+1) - x(ii, jj-1)));
            Ry= ap*(y(ii-1, jj) -2*y(ii, jj) + y(ii+1, jj)) ...
                -0.5*bp*(y(ii+1, jj+1) -y(ii-1, jj+1) - y(ii+1, jj-1) + y(ii-1, jj-1)) ...
                +gp*(y(ii, jj-1) -2*y(ii, jj) + y(ii, jj+1)) + 0.5*dp*(P*(y(ii+1, jj) - y(ii-1,jj)) ...
                + Q*(y(ii, jj+1) - y(ii, jj-1)));
            R = (abs(Rx) + abs(Ry))/2 + R;
        end
    end
    Res = R/((Imax-1)*(Jmax-1));
    if mod(n, 100) == 0
        Res
    end
    ResAr = [ResAr Res];
end

%% 2D EULER FLUID PROBLEM 

U = zeros(4,Imax,Jmax);
F = zeros(4,Imax,Jmax);
G = zeros(4,Imax,Jmax);



