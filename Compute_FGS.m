function [FGS] = Compute_FGS(I_MAX,J_MAX,FG,S)

FGS = zeros(Imax-1,Jmax-1,4,2,4);
for ii = I_MAX-1
    for jj = J_MAX-1
        FGS(ii,jj,1,1,:) = FG(ii,jj,1,:)*S(ii,jj,1,1)/norm(S(ii,jj,1,:));
        FGS(ii,jj,1,2,:) = FG(ii,jj,2,:)*S(ii,jj,1,2)/norm(S(ii,jj,1,:));
        FGS(ii,jj,2,1,:) = FG(ii,jj,1,:)*S(ii,jj,2,1)/norm(S(ii,jj,2,:));
        FGS(ii,jj,2,2,:) = FG(ii,jj,2,:)*S(ii,jj,2,2)/norm(S(ii,jj,2,:));
        FGS(ii,jj,3,1,:) = FG(ii,jj,1,:)*S(ii,jj,3,1)/norm(S(ii,jj,3,:));
        FGS(ii,jj,3,2,:) = FG(ii,jj,2,:)*S(ii,jj,3,2)/norm(S(ii,jj,3,:));
        FGS(ii,jj,4,1,:) = FG(ii,jj,1,:)*S(ii,jj,4,1)/norm(S(ii,jj,4,:));
        FGS(ii,jj,4,2,:) = FG(ii,jj,2,:)*S(ii,jj,4,2)/norm(S(ii,jj,4,:));
    end
end 
end  