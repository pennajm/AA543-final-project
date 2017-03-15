%% AA-543 | JAMES PENNA / CHRISTOPHER PROVENCHER | HW5 | 3/16/2017

clear all; close all; clc

%% NUMERICAL PRESCRIPTION
RES_TOL = 10^(-6); T_TOL = 1; 
I_MAX = 129; J_MAX = 65;
Is = 1:1:I_MAX; Js = 1:1:J_MAX;
x = importdata('xcoors.dat'); y = importdata('ycoors.dat');

%% PHYSICAL PRESCRIPTION
rho_fs = 0.4135; u_fs = 258.4; v_fs = 0.0; C_fs = 304.025; M_fs = 0.85;
E_fs = (1/(1.4-1))*(p_fs/rho_fs) + (u_fs^2 + v_fs^2)/2;

%% SPATIAL DISCRETIZATION
[X, A, S, S_AVG, N, N_AVG] = Spatial_Discretization(I_MAX, J_MAX, x, y);
% X : (I-1,J-1)
% A : (I-1,J-1)
% S : (I-1,J-1,4,2)
% S_AVG : (I-1,J-1,2,2)
% N : (I-1,J-1,4,2)
% N_AVG : (I-1,J-1,2,2)

%% MASTER PROGRAM

U = zeros(I_MAX-1, J_MAX-1,4); 
U(:,:,1) = rho_fs; U(:,:,2) = rho_fs * u_fs; U(:,:,3) = rho_fs * v_fs;
U(:,:,4) = rho_fs * E_fs;

T = 0;
while RES > RES_TOL
    T = T+1;
    if (T > T_TOL)
      disp('BREAK');
      break;
    end
    for RK = 2:4
        FG = Compute_FG(I_MAX,J_MAX,U); % (I-1,J-1,2,4)
        FGS = Compute_FGS(I_MAX,J_MAX,FG,S); % (I-1,J-1,4,2,4)
        FGS_AVG = Compute_FGS_AVG(I_MAX,J_MAX,FGS,S); % (I-1,J-1,4,4)
        D = Compute_D(I_MAX,J_MAX,U,S); % (I-1,J-1,4,4)
        R = Compute_R(I_MAX,J_MAX,D,FGS_AVG); % (I-1,J-1)
        TAU = Compute_TAU(I_MAX,J_MAX,C,U,S_AVG); % (I-1,J-1)
        U = Update_U(I_MAX,J_MAX,RK,U); % (I-1,J-1)
    end
    RES = sum(sum(R));  
end            
            
%% PLOT PHYSICAL QUANTITIES            
  



