function FGS = Give_No_Flux(I_MAX,J_MAX,FG,FGB1,FGB2,S)

FGS = zeros(I_MAX-1,J_MAX-1,4,4);

for ii = 1:(I_MAX-1)
    for jj = 1:(J_MAX-1)     
        if(ii==1)
            FGS(ii,jj,1,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii+1,jj,1,:))*S(ii,jj,1,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii+1,jj,2,:))*S(ii,jj,1,2);
            FGS(ii,jj,2,:) = (1/2)*(FG(ii,jj,1,:) + FG(I_MAX-1,jj,1,:))*S(ii,jj,2,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FG(I_MAX-1,jj,2,:))*S(ii,jj,2,2);               
        elseif(ii==I_MAX-1)
            FGS(ii,jj,1,:) = (1/2)*(FG(ii,jj,1,:) + FG(1,jj,1,:))*S(ii,jj,1,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FG(1,jj,2,:))*S(ii,jj,1,2);
            FGS(ii,jj,2,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii-1,jj,1,:))*S(ii,jj,2,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii-1,jj,2,:))*S(ii,jj,2,2); 
        else
            FGS(ii,jj,1,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii+1,jj,1,:))*S(ii,jj,1,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii+1,jj,2,:))*S(ii,jj,1,2);
            FGS(ii,jj,2,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii-1,jj,1,:))*S(ii,jj,2,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii-1,jj,2,:))*S(ii,jj,2,2);             
        end
        if(jj==1)
            FGS(ii,jj,3,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii,jj+1,1,:))*S(ii,jj,3,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii,jj+1,2,:))*S(ii,jj,3,2);
            FGS(ii,jj,4,:) = (1/2)*(FG(ii,jj,1,:) + FGB1(ii,1,1,:))*S(ii,jj,4,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FGB1(ii,1,2,:))*S(ii,jj,4,2);               
        elseif(jj==J_MAX-1)
            FGS(ii,jj,3,:) = (1/2)*(FG(ii,jj,1,:) + FGB2(ii,1,1,:))*S(ii,jj,3,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FGB2(ii,1,2,:))*S(ii,jj,3,2);
            FGS(ii,jj,4,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii,jj-1,1,:))*S(ii,jj,4,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii,jj-1,2,:))*S(ii,jj,4,2); 
        else
            FGS(ii,jj,3,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii,jj+1,1,:))*S(ii,jj,3,1) ... 
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii,jj+1,2,:))*S(ii,jj,3,2);
            FGS(ii,jj,4,:) = (1/2)*(FG(ii,jj,1,:) + FG(ii,jj-1,1,:))*S(ii,jj,4,1) ...
                           + (1/2)*(FG(ii,jj,2,:) + FG(ii,jj-1,2,:))*S(ii,jj,4,2);             
        end
    end
end
end
  