%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

Imax = 129; Jmax = 65;
Is = 1:1:Imax; Js = 1:1:Jmax;
x = importdata('xcoors.dat'); y = importdata('ycoors.dat');

%% 2D EULER FLUID PROBLEM 

U = zeros(Imax,Jmax,4);
F = zeros(Imax,Jmax,4);
G = zeros(Imax,Jmax,4);

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

%% SPACE DISCRETIZATION

cell_area = zeros(Imax-1,Jmax-1); 
cell_points = zeros(Imax-1,Jmax-1,10); % xy_c, xy_i, xy_i1, xy_j, xy_j1 
cell_v = zeros(Imax-1,Jmax-1,8); % v_i, v_i1, v_j, v_j1
cell_n = zeros(Imax-1,Jmax-1,8); % n_i, n_i1, n_j, n_j1
cell_na = zeros(Imax-1,Jmax-1,4); % n_ai, n_aj
u_IC = zeros(Imax-1,Jmax-1,4); % u_ni, u_nj

for ii = 1:(Imax-1)
    for jj = 1:(Jmax-1)
        % corner vectors
        x_ij1_ij = x(ii,jj+1)-x(ii,jj); y_ij1_ij = y(ii,jj+1)-y(ii,jj);  
        x_i1j1_i1j = x(ii+1,jj+1)-x(ii+1,jj); y_i1j1_i1j = y(ii+1,jj+1)-y(ii+1,jj);
        x_i1j_ij = x(ii+1,jj)-x(ii,jj); y_i1j_ij = y(ii+1,jj)-y(ii,jj);
        x_i1j1_ij1 = x(ii+1,jj+1)-x(ii,jj+1); y_i1j1_ij1 = y(ii+1,jj+1)-y(ii,jj+1);
        % cross vectors
        x_ij1_i1j = x(ii,jj+1)-x(ii+1,jj); y_ij1_i1j = y(ii,jj+1)-y(ii+1,jj); 
        x_i1j1_ij = x(ii+1,jj+1)-x(ii,jj);  y_i1j1_ij = y(ii+1,jj+1)-y(ii,jj); 
        % cell area
        cell_area(ii,jj) = -(1/2)*(x_ij1_i1j * y_i1j1_ij - x_i1j1_ij * y_ij1_i1j);
        % cell points
        cell_points(ii,jj,1) = (1/2)*x_ij1_ij + x(ii,jj);
        cell_points(ii,jj,2) = (1/2)*y_ij1_ij + y(ii,jj);
        cell_points(ii,jj,3) = (1/2)*x_i1j1_i1j + x(ii+1,jj);
        cell_points(ii,jj,4) = (1/2)*y_i1j1_i1j + y(ii+1,jj);
        cell_points(ii,jj,5) = (1/2)*x_i1j_ij + x(ii,jj); 
        cell_points(ii,jj,6) = (1/2)*y_i1j_ij + y(ii,jj);
        cell_points(ii,jj,7) = (1/2)*x_i1j1_ij1 + x(ii,jj+1);
        cell_points(ii,jj,8) = (1/2)*y_i1j1_ij1 + y(ii,jj+1);
        % cell v
        cell_v(ii,jj,1) = x_ij1_ij /((x_ij1_ij^2 + y_ij1_ij^2)^(1/2));
        cell_v(ii,jj,2) = y_ij1_ij /((x_ij1_ij^2 + y_ij1_ij^2)^(1/2));
        cell_v(ii,jj,3) = x_i1j1_i1j /((x_i1j1_i1j^2 + y_i1j1_i1j^2)^(1/2));
        cell_v(ii,jj,4) = y_i1j1_i1j /((x_i1j1_i1j^2 + y_i1j1_i1j^2)^(1/2));
        cell_v(ii,jj,5) = x_i1j_ij /((x_i1j_ij^2 + y_i1j_ij^2)^(1/2));
        cell_v(ii,jj,6) = y_i1j_ij /((x_i1j_ij^2 + y_i1j_ij^2)^(1/2));
        cell_v(ii,jj,7) = x_i1j1_ij1 /((x_i1j1_ij1^2 + y_i1j1_ij1^2)^(1/2));
        cell_v(ii,jj,8) = y_i1j1_ij1 /((x_i1j1_ij1^2 + y_i1j1_ij1^2)^(1/2));
        % cell n  
        cell_n(ii,jj,1) = -cell_v(ii,jj,2); cell_n(ii,jj,2) = cell_v(ii,jj,1);
        cell_n(ii,jj,3) = cell_v(ii,jj,4); cell_n(ii,jj,4) = -cell_v(ii,jj,3);
        cell_n(ii,jj,5) = cell_v(ii,jj,6); cell_n(ii,jj,6) = -cell_v(ii,jj,5);
        cell_n(ii,jj,7) = -cell_v(ii,jj,8); cell_n(ii,jj,8) = cell_v(ii,jj,7);
        % averaged cell n
        cell_na(ii,jj,1) = (1/2)*(cell_n(ii,jj,3) - cell_n(ii,jj,1));
        cell_na(ii,jj,2) = (1/2)*(cell_n(ii,jj,4) - cell_n(ii,jj,2));
        cell_na(ii,jj,3) = (1/2)*(cell_n(ii,jj,7) - cell_n(ii,jj,5));
        cell_na(ii,jj,4) = (1/2)*(cell_n(ii,jj,8) - cell_n(ii,jj,6));
        % u_IC 
        u_IC(ii,jj,1) = (ufs * cell_na(ii,jj,1))*cell_na(ii,jj,1);
        u_IC(ii,jj,2) = (ufs * cell_na(ii,jj,1))*cell_na(ii,jj,2);
        u_IC(ii,jj,3) = (ufs * cell_na(ii,jj,3))*cell_na(ii,jj,3);
        u_IC(ii,jj,4) = (ufs * cell_na(ii,jj,3))*cell_na(ii,jj,4);
        % NOTE: for different angle of attack, change dot product in u_IC code!
    end
end
%pcolor(cell_area)
%quiver(cell_points(:,:,1),cell_points(:,:,2),cell_n(:,:,1),cell_n(:,:,2))
%quiver(cell_points(:,:,3),cell_points(:,:,4),cell_n(:,:,3),cell_n(:,:,4))
%quiver(cell_points(:,:,5),cell_points(:,:,6),cell_n(:,:,5),cell_n(:,:,6))
%quiver(cell_points(:,:,7),cell_points(:,:,8),cell_n(:,:,7),cell_n(:,:,8))
%quiver(cell_points(:,:,1),cell_points(:,:,2),cell_na(:,:,1),cell_na(:,:,2))
%quiver(cell_points(:,:,1),cell_points(:,:,2),cell_na(:,:,3),cell_na(:,:,4))
%quiver(cell_points(:,:,1),cell_points(:,:,2),u_IC(:,:,1),u_IC(:,:,2))





