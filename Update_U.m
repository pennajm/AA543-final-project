function U = Update_U(I_MAX,J_MAX,RK,TAU,R,U)

alpha = zeros(4,1);
alpha(1) = 1/8; alpha(2) = 0.306; alpha(3) = 0.587; alpha(4) = 1;
for ii = 1:I_MAX-1
    for jj = 1:J_MAX-1
        U(ii,jj) = U(ii,jj) - TAU(ii,jj).*alpha(RK).*R(ii,jj);
    end
end
end

