%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

RES_TOL = 10^(-6); T_TOL = 1; CFL = 2.5;
Imax = 129; Jmax = 65;
Is = 1:1:Imax; Js = 1:1:Jmax;
x = importdata('xcoors.dat'); y = importdata('ycoors.dat');

%% FREESTREAM PHYSICAL QUANTITIES 

gamma = 1.4;
rhofs = 0.4135;
ufs = 258.4;
vfs = 0.0;
cfs = 304.025;
Mfs = 0.85;
rhoufs = 0.4135*258.4;
rhovfs = 0.0;
Pfs = 27300.0;
rhoEfs = (1/(1.4-1))*Pfs + 0.5*rhofs*(ufs^2);
rhouH = rhoufs*((gamma/(gamma-1))*Pfs/rhofs + 0.5*rhofs*ufs^2);
rhovH = rhovfs*((gamma/(gamma-1))*Pfs/rhofs + 0.5*rhofs*vfs^2);

Ufs = [rhofs, rhoufs, rhovfs, rhoEfs]';
Ufn = [rhofs, rhoufs, rhoufs, rhoEfs]';
Ufsn = Ufs./[rhofs, rhoufs, rhoufs, rhoEfs]';

Ffs = [rhoufs, (rhofs*ufs^2 +Pfs), rhofs*ufs*vfs, rhouH]';
Gfs = [rhoufs, rhofs*ufs*vfs, (rhofs*vfs^2 + Pfs), rhovH]';

%% SPACE DISCRETIZATION

cell_area = zeros(Imax-1,Jmax-1); 
cell_points = zeros(Imax-1,Jmax-1,10); % xy_c, xy_i, xy_i1, xy_j, xy_j1 
cell_v = zeros(Imax-1,Jmax-1,8); % v_i, v_i1, v_j, v_j1
S = zeros(Imax-1,Jmax-1,8); % n_i, n_i1, n_j, n_j1
S_avg = zeros(Imax-1,Jmax-1,4); % n_ai, n_aj
u = zeros(Imax-1,Jmax-1,4); % u_ni, u_nj
U = zeros(Imax-1, Jmax-1, 4);
F = zeros(Imax-1, Jmax-1, 4,4);
G = zeros(Imax-1, Jmax-1, 4,4);

%% COMPUTE SPATIAL QUANTITIES
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
        % cell n  
        S(ii,jj,1) = -y_ij1_ij; S(ii,jj,2) = x_ij1_ij;
        S(ii,jj,3) = y_i1j1_i1j; S(ii,jj,4) = -x_i1j1_i1j;
        S(ii,jj,5) = y_i1j_ij; S(ii,jj,6) = -x_i1j_ij;
        S(ii,jj,7) = -y_i1j1_ij1; S(ii,jj,8) = x_i1j1_ij1;
        % averaged cell n
        S_avg(ii,jj,1) = (1/2)*(S(ii,jj,3) - S(ii,jj,1));
        S_avg(ii,jj,2) = (1/2)*(S(ii,jj,4) - S(ii,jj,2));
        S_avg(ii,jj,3) = (1/2)*(S(ii,jj,7) - S(ii,jj,5));
        S_avg(ii,jj,4) = (1/2)*(S(ii,jj,8) - S(ii,jj,6));
    end
end

%% COMPUTE PHYSICAL FLUX
for ii = Imax-1
    for jj = Jmax-1

        F(ii,jj,1,:) = (Ffs * S_avg(ii,jj,1))*S_avg(ii,jj,1); %aix
        F(ii,jj,2,:) = (Ffs * S_avg(ii,jj,1))*S_avg(ii,jj,2); %aiy 
        F(ii,jj,3,:) = (Ffs * S_avg(ii,jj,3))*S_avg(ii,jj,3); %ajx
        F(ii,jj,4,:) = (Ffs * S_avg(ii,jj,3))*S_avg(ii,jj,4); %ajy
        
        G(ii,jj,1,:) = (Gfs * S_avg(ii,jj,2))*S_avg(ii,jj,1); %aix
        G(ii,jj,2,:) = (Gfs * S_avg(ii,jj,2))*S_avg(ii,jj,2); %aiy 
        G(ii,jj,3,:) = (Gfs * S_avg(ii,jj,4))*S_avg(ii,jj,3); %ajx
        G(ii,jj,4,:) = (Gfs * S_avg(ii,jj,4))*S_avg(ii,jj,4); %ajy
        % NOTE: for different angle of attack, change dot product in u_IC code!
    end
end
%% COMPUTE ARTIFICIAL VISCOSITY
F_flux = zeros(Imax-1, Jmax-1, 4,4);
G_flux = zeros(Imax-1, Jmax-1, 4,4);
R = zeros(Imax-1, Jmax-1);
D = zeros(Imax-1, Jmax-1, 4, 4); 

Rnbp = u(round(0.25*Imax):round(0.75*Imax),Jmax-1,3) + ...
     + u(round(0.25*Imax):round(0.75*Imax), Jmax-1,4) + 2*cfs/(gamma-1);
Rnim = u(round(0.25*Imax):round(0.75*Imax),Jmax-1,3) + ...
     +u(round(0.25*Imax):round(0.75*Imax), Jmax-1,4) - 2*cfs/(gamma-1);
%unbi = (Rnbp + Rnim)/2;
cbi = (Rnbp - Rnim)*0.25*(gamma-1);

P = (gamma -1)*(U(:,:,4) - 0.5*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
uvec = ones(Imax,Jmax,2);
uvec(:,:,1) = ufs; 
uvec(:,:,2) = vfs;
nu = zeros(Imax,1);
E2 = zeros(Imax-1,Jmax-1,4);
E4 = zeros(Imax-1,Jmax-1,4);
k2 = 1/4; 
k4 = 1/256;

for ii = 1:Imax-1
 
    if ii==1
        nu(ii,1) = abs((P(ii+1,:) - 2*P(ii,:) + P(Imax-1,:))/...
                (P(ii+1,:) + 2*P(ii,:) + P(Imax-1,:)));
    elseif ii ==Imax-1
        nu(ii,1) = abs((P(1,:) - 2*P(ii,:) + P(ii-1,:))/...
        (P(1,:) + 2*P(ii,:) + P(ii-1,:)));
    else
        nu(ii,1) = abs((P(ii+1,:) - 2*P(ii,:) + P(ii-1,:))/...
        (P(ii+1,:) + 2*P(ii,:) + P(ii-1,:)));
    end
end

for nind = 0:3
    for ii = 1:Imax-1
        for jj = 1:Jmax-1
            magds = sqrt(cell_n(ii,jj,(1+nind*2))^2 + cell_n(ii,jj,(2+nind*2))^2);
            if ii == 1
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + cfs*magds)*max([nu(Imax),nu(ii), nu(ii+1), nu(ii+2)]);
            elseif ii == Imax-1
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + cfs*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(1)]);
            else
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + cfs*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(ii+2)]);
            end

            E4(ii,jj,nind+1) = max([0,(0.5*k4*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                + cfs*magds) - E2(ii,jj,nind+1))]);

            if ii== 1
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(ii+2,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:)...
                    + 3*U(ii,jj,nind+1,:) - U(Imax-1,jj,nind+1,:)));
            elseif ii == Imax-2
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(1,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:) + ...
                    3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));

            elseif ii == Imax-1
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(2,jj,nind+1,:) - 3*U(1,jj,nind+1,:) + ...
                    3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));
            else
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj, nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(ii+2,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:)...
                    + 3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));
            end
        end
    end
end

%% COMPUTE NUMERICAL FLUX / COMPUTE RESIDUAL

for ii = 1:(Imax-1)
    for jj = 1:(Jmax-1)     
        if ii < 128
            F_flux(ii,jj,1,:) = (1/2)*(F(ii,jj,1,:)+F(ii+1,jj,1,:))*cell_n(ii,jj,1)-D(ii,jj,1,:);
            F_flux(ii,jj,2,:) = (1/2)*(F(ii,jj,2,:)+F(ii+1,jj,2,:))*cell_n(ii,jj,3)-D(ii,jj,2,:);
            F_flux(ii,jj,3,:) = (1/2)*(F(ii,jj,3,:)+F(ii+1,jj,3,:))*cell_n(ii,jj,5)-D(ii,jj,3,:);
            F_flux(ii,jj,4,:) = (1/2)*(F(ii,jj,4,:)+F(ii+1,jj,4,:))*cell_n(ii,jj,7)-D(ii,jj,4,:);
            G_flux(ii,jj,1,:) = (1/2)*(G(ii,jj,1,:)+G(ii+1,jj,1,:))*cell_n(ii,jj,2)-D(ii,jj,1,:);
            G_flux(ii,jj,2,:) = (1/2)*(G(ii,jj,2,:)+G(ii+1,jj,2,:))*cell_n(ii,jj,4)-D(ii,jj,2,:);
            G_flux(ii,jj,3,:) = (1/2)*(G(ii,jj,3,:)+G(ii+1,jj,3,:))*cell_n(ii,jj,5)-D(ii,jj,3,:);
            G_flux(ii,jj,4,:) = (1/2)*(G(ii,jj,4,:)+G(ii+1,jj,4,:))*cell_n(ii,jj,8)-D(ii,jj,4,:);
            sum_F_flux = sum(F_flux(ii,jj,1,:)+F_flux(ii,jj,2,:)+F_flux(ii,jj,3,:)+F_flux(ii,jj,4,:));
            sum_G_flux = sum(G_flux(ii,jj,1,:)+G_flux(ii,jj,2,:)+G_flux(ii,jj,3,:)+G_flux(ii,jj,4,:));
            R(ii,jj) = sum_F_flux + sum_G_flux;
        else
            F_flux(ii,jj,1,:) = (1/2)*(F(ii,jj,1,:)+F(1,jj,1,:))*cell_n(ii,jj,1)-D(ii,jj,1,:);
            F_flux(ii,jj,2,:) = (1/2)*(F(ii,jj,2,:)+F(1,jj,2,:))*cell_n(ii,jj,3)-D(ii,jj,2,:);
            F_flux(ii,jj,3,:) = (1/2)*(F(ii,jj,3,:)+F(1,jj,3,:))*cell_n(ii,jj,5)-D(ii,jj,3,:);
            F_flux(ii,jj,4,:) = (1/2)*(F(ii,jj,4,:)+F(1,jj,4,:))*cell_n(ii,jj,7)-D(ii,jj,4,:);
            G_flux(ii,jj,1,:) = (1/2)*(G(ii,jj,1,:)+G(1,jj,1,:))*cell_n(ii,jj,2)-D(ii,jj,1,:);
            G_flux(ii,jj,2,:) = (1/2)*(G(ii,jj,2,:)+G(1,jj,2,:))*cell_n(ii,jj,4)-D(ii,jj,2,:);
            G_flux(ii,jj,3,:) = (1/2)*(G(ii,jj,3,:)+G(1,jj,3,:))*cell_n(ii,jj,5)-D(ii,jj,3,:);
            G_flux(ii,jj,4,:) = (1/2)*(G(ii,jj,4,:)+G(1,jj,4,:))*cell_n(ii,jj,8)-D(ii,jj,4,:);
            sum_F_flux = sum(F_flux(ii,jj,1,:)+F_flux(ii,jj,2,:)+F_flux(ii,jj,3,:)+F_flux(ii,jj,4,:));
            sum_G_flux = sum(G_flux(ii,jj,1,:)+G_flux(ii,jj,2,:)+G_flux(ii,jj,3,:)+G_flux(ii,jj,4,:));
            R(ii,jj) = sum_F_flux + sum_G_flux;
        end
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

% %from here assuming the normal boundary conditions
% U = U_IC; G = G_IC; F = F_IC;
% %calculate Reimann invariants based on freestream
% Rnbpi = u_IC(0.25*Imax:0.75*Imax,Jmax-1,3) + ...
%     +u_IC(0.25*Imax:0.75*Imax, Jmax-1,4) + 2*cfs/(gamma-1);
% Rnim = u_IC(0.25*Imax:0.75*Imax,Jmax-1,3) + ...
%     +u_IC(0.25*Imax:0.75*Imax, Jmax-1,4) - 2*cfs/(gamma-1);
% %unbi = (Rnbp + Rnim)/2;
% cbi = (Rnbp - Rnim)*0.25*(gamma-1);
% %use streamline constants to calculate interior values
% Hfs = (gamma/(gamma-1))*Pfs/rhofs + 0.5*ufs^2;


%% COMPUTE TIME STEPS

alpha = zeros(4,1);
alpha(1) = 1/8; alpha(2) = 0.306; alpha(3) = 0.587; alpha(4) = 1;
tau = zeros(ii,jj);
for ii = 1:Imax-1
    for jj = 1:Jmax-1
        shock_i = (U(ii,jj,2)/U(ii,jj,1)+C(ii,jj))*S_avg(ii,jj,1) ...
                + (U(ii,jj,3)/U(ii,jj,1)+C(ii,jj))*S_avg(ii,jj,2);
        shock_j = (U(ii,jj,2)/U(ii,jj,1)+C(ii,jj))*S_avg(ii,jj,3) ...
                + (U(ii,jj,3)/U(ii,jj,1)+C(ii,jj))*S_avg(ii,jj,4);
        tau(ii,jj) = CFL/(abs(shock_i)+abs(shock_j));
    end
end

%% TIME EVOLUTION
%T = 0;
% while RES < RES_TOL
%     T = T+1;
%     if (T > T_TOL)
%       disp('BREAK');
%       break;
%     end
%     % COMPUTE TAU
%     for RK = 2:4
%         % COMPUTE F
%         % COMPUTE D
%         % COMPUTE F_FLUX, G_FLUX, RESIDUAL
%         for ii = 1:Imax-1
%             for jj = 1:Jmax-1
%                 U(ii,jj) = U(ii,jj) - tau(ii,jj)*alpha(RK)*R(ii,jj);
%             end
%         end
%     end
%     RES = sum(sum(R));
%      
% end            
            
            
  



