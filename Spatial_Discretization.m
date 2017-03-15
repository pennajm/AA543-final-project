function [X, A, S, S_AVG, N, N_AVG] = Spatial_Discretization(I_MAX, J_MAX, x, y)

X = zeros(I_MAX-1,J_MAX-1);
A = zeros(I_MAX-1,J_MAX-1);
S = zeros(I_MAX-1,J_MAX-1,4,2);
S_AVG = zeros(I_MAX-1,J_MAX-1,2,2);
N = zeros(I_MAX-1,J_MAX-1,4,2);
N_AVG = zeros(I_MAX-1,J_MAX-1,2,2);

for ii = 1:(I_MAX-1)
    for jj = 1:(J_MAX-1)
        % corner vectors
        x_ij1_ij = x(ii,jj+1)-x(ii,jj); y_ij1_ij = y(ii,jj+1)-y(ii,jj);  
        x_i1j1_i1j = x(ii+1,jj+1)-x(ii+1,jj); y_i1j1_i1j = y(ii+1,jj+1)-y(ii+1,jj);
        x_i1j_ij = x(ii+1,jj)-x(ii,jj); y_i1j_ij = y(ii+1,jj)-y(ii,jj);
        x_i1j1_ij1 = x(ii+1,jj+1)-x(ii,jj+1); y_i1j1_ij1 = y(ii+1,jj+1)-y(ii,jj+1);
        % cross vectors
        x_ij1_i1j = x(ii,jj+1)-x(ii+1,jj); y_ij1_i1j = y(ii,jj+1)-y(ii+1,jj); 
        x_i1j1_ij = x(ii+1,jj+1)-x(ii,jj);  y_i1j1_ij = y(ii+1,jj+1)-y(ii,jj); 
        
        A(ii,jj) = -(1/2)*(x_ij1_i1j * y_i1j1_ij - x_i1j1_ij * y_ij1_i1j); 
        
        S(ii,jj,1,1) = -y_ij1_ij; S(ii,jj,1,2) = x_ij1_ij;
        S(ii,jj,2,1) = y_i1j1_i1j; S(ii,jj,2,2) = -x_i1j1_i1j;
        S(ii,jj,3,1) = y_i1j_ij; S(ii,jj,3,2) = -x_i1j_ij;
        S(ii,jj,4,1) = -y_i1j1_ij1; S(ii,jj,4,2) = x_i1j1_ij1;
        
        S_AVG(ii,jj,1,1) = (1/2)*(S(ii,jj,2,1) - S(ii,jj,1,1));
        S_AVG(ii,jj,1,2) = (1/2)*(S(ii,jj,2,2) - S(ii,jj,1,2));
        S_AVG(ii,jj,2,1) = (1/2)*(S(ii,jj,4,1) - S(ii,jj,3,1));
        S_AVG(ii,jj,2,2) = (1/2)*(S(ii,jj,4,2) - S(ii,jj,3,2));
        
        N(ii,jj,1,:) = S(ii,jj,1,:)/norm(S(ii,jj,1,:));
        N(ii,jj,2,:) = S(ii,jj,2,:)/norm(S(ii,jj,2,:));
        N(ii,jj,3,:) = S(ii,jj,3,:)/norm(S(ii,jj,3,:));
        N(ii,jj,4,:) = S(ii,jj,4,:)/norm(S(ii,jj,4,:));
        
        N_AVG(ii,jj,1,:) = S_AVG(ii,jj,1,:)/norm(S_AVG(ii,jj,1,:));
        N_AVG(ii,jj,2,:) = S_AVG(ii,jj,2,:)/norm(S_AVG(ii,jj,2,:));
        N_AVG(ii,jj,3,:) = S_AVG(ii,jj,3,:)/norm(S_AVG(ii,jj,3,:));
        N_AVG(ii,jj,4,:) = S_AVG(ii,jj,4,:)/norm(S_AVG(ii,jj,4,:));
        
    end
end


end

