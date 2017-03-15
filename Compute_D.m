function D = Compute_D(I_MAX,J_MAX,U,S)

P = (1.4-1)*(U(:,:,4) - (1/2)*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
C = (1.4*P./U(:,:,1));

u = U(:,:,2)./U(:,:,1); 
v = U(:,:,3)./U(:,:,1);
nu = zeros(I_MAX-1,1);
E2 = zeros(I_MAX-1,J_MAX-1,4);
E4 = zeros(I_MAX-1,J_MAX-1,4);
k2 = 1/4; 
k4 = 1/256;
D = zeros(I_MAX-1,J_MAX-1,4);

for ii = 1:I_MAX-1
    if ii==1
        nu(ii,1) = abs((P(ii+1,:) - 2*P(ii,:) + P(I_MAX-1,:))/...
                (P(ii+1,:) + 2*P(ii,:) + P(I_MAX-1,:)));
    elseif ii ==I_MAX-1
        nu(ii,1) = abs((P(1,:) - 2*P(ii,:) + P(ii-1,:))/...
        (P(1,:) + 2*P(ii,:) + P(ii-1,:)));
    else
        nu(ii,1) = abs((P(ii+1,:) - 2*P(ii,:) + P(ii-1,:))/...
        (P(ii+1,:) + 2*P(ii,:) + P(ii-1,:)));
    end
end

for nind = 1:4
    for ii = 1:I_MAX-1
        for jj = 1:J_MAX-1
            magds = sqrt(S(ii,jj,nind,1)^2 + S(ii,jj,nind,2)^2);
            if ii == 1
                E2(ii,jj,nind) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C*magds)*max([nu(I_MAX-1),nu(ii), nu(ii+1), nu(ii+2)]);
            elseif ii == I_MAX-1
                E2(ii,jj,nind) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(1)]);
            else
                E2(ii,jj,nind) = 0.5*k2*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                    + C*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(ii+2)]);
            end

            E4(ii,jj,nind) = max([0,(0.5*k4*((u(ii,jj)*S(ii,jj,nind,1) + v(ii,jj)*S(ii,jj,nind,2))...
                + C*magds) - E2(ii,jj,nind))]);

            if ii== 1
                D(ii,jj,nind) = E2(ii,jj,nind).*(U(ii+1,jj) - U(ii,jj)) - ...
                    (E4(ii,jj,nind).*(U(ii+2,jj) - 3*U(ii+1,jj)...
                    + 3*U(ii,jj) - U(I_MAX-1,jj)));
            elseif ii == I_MAX-2
                D(ii,jj,nind) = E2(ii,jj,nind).*(U(ii+1,jj) - U(ii,jj)) - ...
                    (E4(ii,jj,nind).*(U(1,jj) - 3*U(ii+1,jj) + ...
                    3*U(ii,jj) - U(ii-1,jj)));

            elseif ii == I_MAX-1
                D(ii,jj,nind) = E2(ii,jj,nind).*(U(1,jj) - U(ii,jj)) - ...
                    (E4(ii,jj,nind).*(U(2,jj) - 3*U(1,jj) + ...
                    3*U(ii,jj) - U(ii-1,jj)));
            else
                D(ii,jj,nind) = E2(ii,jj,nind).*(U(ii+1,jj) - U(ii,jj)) - ...
                    (E4(ii,jj,nind).*(U(ii+2,jj) - 3*U(ii+1,jj)...
                    + 3*U(ii,jj) - U(ii-1,jj)));
            end
        end
    end
end


end

