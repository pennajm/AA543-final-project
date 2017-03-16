function TAU = Compute_TAU(I_MAX,J_MAX,U,S_AVG)

CFL = 1.0;
TAU = zeros(I_MAX-1,J_MAX-1);
for ii = 1:I_MAX-1
    for jj = 1:J_MAX-1
        C = (1.4 * ((1.4-1)/U(ii,jj,1))*(U(ii,jj,4)...
          -(1/2)*((U(ii,jj,2)^2 + U(ii,jj,3)^2)/U(ii,jj,1))))^(1/2);
        TAU(ii,jj) = CFL/(abs((U(ii,jj,2)/U(ii,jj,1)+C)*S_AVG(ii,jj,1,1) ...
                   + (U(ii,jj,3)/U(ii,jj,1)+C)*S_AVG(ii,jj,1,2)) ...
                   + abs((U(ii,jj,2)/U(ii,jj,1)+C)*S_AVG(ii,jj,2,1) ...
                   + (U(ii,jj,3)/U(ii,jj,1)+C)*S_AVG(ii,jj,2,2)));
    end
end
%TAU = TAU ./ A;
end

