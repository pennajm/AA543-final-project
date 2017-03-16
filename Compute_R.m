function R = Compute_R(I_MAX,J_MAX,D,FGS)

R = zeros(I_MAX-1,J_MAX-1,4);
    for ii = 1:I_MAX-1
        for jj=1:J_MAX-1
            R(ii,jj,1) = sum(FGS(ii,jj,:,1))-D(ii,jj,1);
            R(ii,jj,2) = sum(FGS(ii,jj,:,2))-D(ii,jj,2);
            R(ii,jj,3) = sum(FGS(ii,jj,:,3))-D(ii,jj,3);
            R(ii,jj,4) = sum(FGS(ii,jj,:,4))-D(ii,jj,4);
        end
    end
end

