%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

Imax = 129; Jmax = 65;
Is = 1:1:Imax; Js = 1:1:Jmax;
x = importdata('grid_x.dat'); y = importdata('grid_y.dat');


%% 2D EULER FLUID PROBLEM 

U = zeros(Imax,Jmax,4);
F = zeros(Imax,Jmax,4);
G = zeros(Imax,Jmax,4);
cells = zeros(Imax-1,Jmax-1,5); 

%% SPACE DISCRETIZATION
for ii = 1:(Imax-1)
    for jj = 1:(Jmax-1)
        cells(ii,jj,1) = (1/4)*(x(ii,jj)+x(ii+1,jj)+x(ii,jj)+x(ii,jj+1)); % x_center
        cells(ii,jj,2) = (1/4)*(y(ii,jj)+y(ii+1,jj)+y(ii,jj)+y(ii,jj+1)); % y_center
        x_ac = x(ii+1,jj) - x(ii,jj+1); y_ac = y(ii+1,jj) - y(ii,jj+1);
        x_bd = x(ii+1,jj+1) - x(ii,jj); y_bd = y(ii+1,jj+1) - y(ii,jj);
        cells(ii,jj,3) = (1/2)*(x_ac * y_bd - x_bd * y_ac); % area
        x_ba = x(ii+1,jj+1)-x(ii+1,jj); y_ba = y(ii+1,jj+1)-y(ii+1,jj);
        m_ba = (x_ba^2 + y_ba^2)^(1/2); % magnitude of (x_ba,y_ba)
        l_x = x_ba / m_ba; l_y = y_ba / m_ba; % unit vector
        s_x = y_ba*l_x; s_y = -x_ba * l_y; % normal vector
        cells(ii,jj,4) = s_x;
        cells(ii,jj,5) = s_y;
    end
end
%% FREESTREAM PHYSICAL QUANTITIES 

rhofs = 0.4135;
ufs = 258.4;
cfs = 304.025;
Mfs = 0.85;
rhoufs = 0.4135*258.4;
rhovfs = 0.0;
Pfs = 27300.0;
rhoEfs = (1/(1.4-1))*Pfs + 0.5*rhofs*(ufs^2);

Ufs = [rhofs, rhoufs, rhovfs, rhoEfs];
Ufn = [rhofs, rhoufs, rhoufs, rhoEfs];
Ufsn = Ufs./[rhofs, rhoufs, rhoufs, rhoEfs];




