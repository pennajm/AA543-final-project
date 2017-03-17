%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

%% NUMERICAL PRESCRIPTION
RES_TOL = 10^(-6); T_TOL = 5; 
I_MAX = 129; J_MAX = 65;
x = importdata('xcoors.dat'); y = importdata('ycoors.dat');

%% PHYSICAL PRESCRIPTION
rho_fs = 0.4135; u_fs = 258.4; v_fs = 0.0; C_fs = 304.025; M_fs = 0.85;
p_fs = 27300.0;
E_fs = (1/(1.4-1))*(p_fs/rho_fs) + (u_fs^2 + v_fs^2)/2;

%% NUMERICAL DISCRETIZATION
[X, A, S, S_AVG] = Spatial_Discretization(I_MAX, J_MAX, x, y);

%% PHYSICAL INITIALIZATION

U = zeros(I_MAX-1, J_MAX-1,4); 
U(:,:,1) = rho_fs; U(:,:,2) = rho_fs * u_fs; U(:,:,3) = rho_fs * v_fs;
U(:,:,4) = rho_fs * E_fs;

%% MAIN CODE

T = 0; RES = 1; RES_FS = 1;
while (RES/RES_FS) > RES_TOL
    T = T+1;
    if (T > T_TOL)
      disp('BREAK');
      break;
    end
    U_RK = U;
    UB1_RK = Compute_UB1(I_MAX, J_MAX, U,S,X);
    UB2_RK = Compute_UB2(I_MAX, J_MAX, U,S);
    TAU = Compute_TAU(I_MAX,J_MAX,U,S_AVG);
    for RK = 1:4
            [FG,FGB1,FGB2] = Compute_FG(I_MAX,J_MAX,U_RK,UB1_RK,UB2_RK);
            FGS = Give_No_Flux(I_MAX,J_MAX,FG,FGB1,FGB2,S);
            D = Compute_D(I_MAX,J_MAX,U_RK,S,UB1_RK,UB2_RK);    
            R = Compute_R(I_MAX,J_MAX,D,FGS); 
            U_RK = Update_U(I_MAX,J_MAX,RK,TAU,R,U); 
    end
    if(T == 1)
        RES_FS = abs(mean(mean(mean(R))));
    end
    RES = abs(mean(mean(mean(R))));
    RES/RES_FS
    U = U_RK;
end            

%% PLOT PHYSICAL QUANTITIES            

u_mag = ((U(:,:,2)./U(:,:,1)).^2 + (U(:,:,3)./U(:,:,1)).^2).^(1/2);
p = (1.4-1)*(U(:,:,4)-(1/2)*((U(:,:,2).^2 + U(:,:,3).^2)./(U(:,:,1).^2)));
H = U(:,:,4)./U(:,:,1);

% figure
% contour(X(:,:,1),X(:,:,2),u_mag)
% title('Velocity Magnitude')
% xlabel('X')
% ylabel('Y')
% figure
% contour(X(:,:,1),X(:,:,2),p)
% title('Pressure')
% xlabel('X')
% ylabel('Y')
% figure
% contour(X(:,:,1),X(:,:,2),H)
% title('Enthalpy')
% xlabel('X')
% ylabel('Y')



