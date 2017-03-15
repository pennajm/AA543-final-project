%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

%% NUMERICAL PRESCRIPTION
RES_TOL = 10^(-6); T_TOL = 1; 
I_MAX = 129; J_MAX = 65;
x = importdata('xcoors.dat'); y = importdata('ycoors.dat');

%% PHYSICAL PRESCRIPTION
rho_fs = 0.4135; u_fs = 258.4; v_fs = 0.0; C_fs = 304.025; M_fs = 0.85;
p_fs = 27300.0;
E_fs = (1/(1.4-1))*(p_fs/rho_fs) + (u_fs^2 + v_fs^2)/2;

%% SPATIAL DISCRETIZATION
[X, A, S, S_AVG] = Spatial_Discretization(I_MAX, J_MAX, x, y);

%% INITIAL CONDITIONS

U = zeros(I_MAX-1, J_MAX-1,4); 
U(:,:,1) = rho_fs; U(:,:,2) = rho_fs * u_fs; U(:,:,3) = rho_fs * v_fs;
U(:,:,4) = rho_fs * E_fs;

%% BOUNDARY CONDITIONS

T = 0; RES = 1000;
while RES > RES_TOL
    T = T+1;
    if (T > T_TOL)
      disp('BREAK');
      break;
    end
    for RK = 1:4
        UB1 = Compute_UB1(I_MAX, J_MAX, U,S,X);
        UB2 = Compute_UB2(I_MAX, J_MAX, U,S);
        [FG,FGB1,FGB2] = Compute_FG(I_MAX,J_MAX,U,UB1,UB2); 
        [FGS,FGSB1,FGSB2] = Compute_FGS(I_MAX,J_MAX,FG,FGB1,FGB2,S); 
        FGS_AVG = Compute_FGS_AVG(I_MAX,J_MAX,FGS,FGSB1,FGSB2,S); 
        D = Compute_D(I_MAX,J_MAX,U,S);
        R = Compute_R(I_MAX,J_MAX,D,FGS_AVG); 
        TAU = Compute_TAU(I_MAX,J_MAX,U,S_AVG); 
        U = Update_U(I_MAX,J_MAX,RK,TAU,R,U); 
    end
    %RES = sum(sum(R));  
end            
            
%% PLOT PHYSICAL QUANTITIES            
  


