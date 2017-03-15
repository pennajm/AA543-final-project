function FGS_AVG  = Compute_FGS_AVG(I_MAX,J_MAX,FGS,S)

FGS_AVG = zeros(Imax-1,Jmax-1,4,4);
for ii = 1:(I_MAX-1)
    for jj = 1:(J_MAX-1)     
        if(ii < 128) 
            FGS_AVG(ii,jj,1,:) = (1/2)*(FGS(ii,jj,1,1,:)+FGS(ii+1,jj,1,1,:))*S(ii,jj,1,1) ... 
                               + (1/2)*(FGS(ii,jj,1,2,:)+FGS(ii+1,jj,1,2,:))*S(ii,jj,1,2);
        else
            FGS_AVG(ii,jj,1,:) = (1/2)*(FGS(ii,jj,1,1,:)+FGS(ii+1,jj,1,1,:))*S(ii,jj,1,1) ... 
                               + (1/2)*(FGS(ii,jj,1,2,:)+FGS(ii+1,jj,1,2,:))*S(ii,jj,1,2);
        end    
        FGS_AVG(ii,jj,2,:) = (1/2)*(FGS(ii,jj,2,1,:)+FGS(ii-1,jj,2,1,:))*S(ii,jj,1,1) ...
                           + (1/2)*(FGS(ii,jj,2,2,:)+FGS(ii-1,jj,2,2,:))*S(ii,jj,1,2);
        FGS_AVG(ii,jj,3,:) = (1/2)*(FGS(ii,jj,3,1,:)+FGS(ii-1,jj,3,1,:))*S(ii,jj,1,1) ...
                           + (1/2)*(FGS(ii,jj,3,2,:)+FGS(ii-1,jj,3,2,:))*S(ii,jj,1,2); 
        FGS_AVG(ii,jj,4,:) = (1/2)*(FGS(ii,jj,4,1,:)+FGS(ii-1,jj,4,1,:))*S(ii,jj,1,1) ...
                           + (1/2)*(FGS(ii,jj,4,2,:)+FGS(ii-1,jj,4,2,:))*S(ii,jj,1,2);                          
    end
end
end

