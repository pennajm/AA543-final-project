function D = Compute_D(I_MAX,J_MAX,U,S, UB1, UB2)

P = (1.4-1)*(U(:,:,4) - (1/2)*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
PB1 = (1.4-1)*(UB1(:,4) - (1/2)*(UB1(:,2).^2 + UB1(:,3).^2)./UB1(:,1));
PB2 = (1.4-1)*(UB2(:,4) - (1/2)*(UB2(:,2).^2 + UB2(:,3).^2)./UB2(:,1));
C = (1.4*P./U(:,:,1));

u = U(:,:,2)./U(:,:,1);
v = U(:,:,3)./U(:,:,1);
nui = zeros(I_MAX-1,J_MAX-1);
nuj = zeros(I_MAX-1,J_MAX-1);
E2 = zeros(I_MAX-1,J_MAX-1,4,4);
E4 = zeros(I_MAX-1,J_MAX-1,4,4);
k2 = 1/4;
k4 = 1/256;
D = zeros(I_MAX-1,J_MAX-1,4,4);

for ii = 1:I_MAX-1
    for jj = 1:J_MAX-1
        if ii==1
            nui(ii,jj) = abs((P(ii+1,jj) - 2*P(ii,jj) + P(I_MAX-1,jj))/...
                (P(ii+1,jj) + 2*P(ii,jj) + P(I_MAX-1,jj)));
        elseif ii ==I_MAX-1
            nui(ii,jj) = abs((P(1,jj) - 2*P(ii,jj) + P(ii-1,jj))/...
                (P(1,jj) + 2*P(ii,jj) + P(ii-1,jj)));
        else
            nui(ii,jj) = abs((P(ii+1,jj) - 2*P(ii,jj) + P(ii-1,jj))/...
                (P(ii+1,jj) + 2*P(ii,jj) + P(ii-1,jj)));
        end
    end
end

for jj = 1:J_MAX-1
    for ii = 1:I_MAX-1
        if jj==1
            nuj(ii,jj) = abs((P(ii,jj+1) - 2*P(ii,jj) + PB1(ii))/...
                (P(ii,jj+1) + 2*P(ii,jj) + PB1(ii)));
        elseif jj ==J_MAX-1
            nuj(ii,jj) = abs((PB2(ii) - 2*P(ii,jj) + P(ii,jj-1))/...
                (PB2(ii) + 2*P(ii,jj) + P(ii,jj-1)));
        else
            nuj(ii,jj) = abs((P(ii,jj+1) - 2*P(ii,jj) + P(ii,jj-1))/...
                (P(ii,jj+1) + 2*P(ii,jj) + P(ii,jj-1)));
        end
    end
end


for nind = 1:2
    for ii = 1:I_MAX-1
        for jj = 1:J_MAX-1
            magds = sqrt(S(ii,jj,nind,1)^2 + S(ii,jj,nind,2)^2);
            if ii == 1
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nui(I_MAX-1), nui(ii), nui(ii+1), nui(ii+2)]);
            elseif ii == I_MAX-2
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nui(ii-1),nui(ii), nui(ii+1), nui(1)]);
            elseif ii == I_MAX-1
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nui(ii-1),nui(ii), nui(1), nui(2)]);
            else
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nui(ii-1),nui(ii), nui(ii+1), nui(ii+2)]);
            end
            
        end
    end
end
for nind = 3:4
    for ii = 1:I_MAX-1
        for jj = 1:J_MAX-1
            magds = sqrt(S(ii,jj,nind,1)^2 + S(ii,jj,nind,2)^2);
            if jj == 1
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nuj(I_MAX-1),nuj(jj), nuj(jj+1), nuj(jj+2)]);
            elseif jj == J_MAX-2
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nuj(jj-1),nuj(jj), nuj(jj+1), nuj(J_MAX-1)]);
            elseif jj == J_MAX-1
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nuj(jj-1),nuj(jj), nuj(J_MAX-1), nuj(J_MAX-1)]);
            else
                E2(ii,jj,nind,:) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C(ii,jj)*magds)*max([nuj(jj-1),nuj(jj), nuj(jj+1), nuj(jj+2)]);
            end
            
        end
    end
end
for nind =1:4
for ii=1:I_MAX-1
    for jj=1:J_MAX-1
        fac = 0.5*k4*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
            + C(ii,jj)*magds);
        E4(ii,jj,nind,1) = max([0,(-E2(ii,jj,nind,1)+fac)]);
        E4(ii,jj,nind,2) = max([0,(-E2(ii,jj,nind,2)+fac)]);
        E4(ii,jj,nind,3) = max([0,(-E2(ii,jj,nind,3)+fac)]);
        E4(ii,jj,nind,4) = max([0,(-E2(ii,jj,nind,4)+fac)]);
        for comp =1:4
            if ii== 1
                D(ii,jj,nind,comp) = E2(ii,jj,nind,comp).*(U(ii+1,jj,comp) - U(ii,jj,comp)) - ...
                    (E4(ii,jj,nind,comp).*(U(ii+2,jj,comp) - 3*U(ii+1,jj,comp)...
                    + 3*U(ii,jj,comp) - U(I_MAX-1,jj,comp)));
            elseif ii == I_MAX-2
                D(ii,jj,nind,comp) = E2(ii,jj,nind,comp).*(U(ii+1,jj,comp) - U(ii,jj,comp)) - ...
                    (E4(ii,jj,nind,comp).*(U(1,jj,comp) - 3*U(ii+1,jj,comp) + ...
                    3*U(ii,jj,comp) - U(ii-1,jj,comp)));
                
            elseif ii == I_MAX-1
                D(ii,jj,nind,comp) = E2(ii,jj,nind,comp).*(U(1,jj,comp) - U(ii,jj,comp)) - ...
                    (E4(ii,jj,nind,comp).*(U(2,jj,comp) - 3*U(1,jj,comp) + ...
                    3*U(ii,jj,comp) - U(ii-1,jj,comp)));
            else
                D(ii,jj,nind,comp) = E2(ii,jj,nind,comp).*(U(ii+1,jj,comp) - U(ii,jj,comp)) - ...
                    (E4(ii,jj,nind,comp).*(U(ii+2,jj,comp) - 3*U(ii+1,jj,comp)...
                    + 3*U(ii,jj,comp) - U(ii-1,jj,comp)));
            end
        end
    end
end
end

D = abs(D);

end

