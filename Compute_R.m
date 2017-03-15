function R = Compute_R(I_MAX,J_MAX,D,FGS_AVG)

R = zeros(I_MAX-1,J_MAX-1);
    for ii = 1:I_MAX-1
        for jj=1:J_MAX-1
            R(ii,jj) = sum(sum(FGS_AVG(ii,jj,:,:)-D(ii,jj,:,:)));
        end
    end
end

