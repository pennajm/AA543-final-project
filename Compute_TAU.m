function TAU = Compute_TAU(I_MAX,J_MAX,C,U,S_AVG)

CFL = 2.5;
TAU = zeros(I_MAX-1,J_MAX-1);
for ii = 1:I_MAX-1
    for jj = 1:J_MAX-1
        TAU(ii,jj) = CFL/(abs((U(ii,jj,2)/U(ii,jj,1)+C(ii,jj))*S_AVG(ii,jj,1) ...
                   + (U(ii,jj,3)/U(ii,jj,1)+C(ii,jj))*S_AVG(ii,jj,2)) ...
                   + abs((U(ii,jj,2)/U(ii,jj,1)+C(ii,jj))*S_AVG(ii,jj,3) ...
                   + (U(ii,jj,3)/U(ii,jj,1)+C(ii,jj))*S_AVG(ii,jj,4)));
    end
end
end

