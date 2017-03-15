function [ output_args ] = Untitled13( input_args )

Rnbp = u(round(0.25*I_MAX):round(0.75*I_MAX),J_MAX-1,3) + ...
     + u(round(0.25*I_MAX):round(0.75*I_MAX), J_MAX-1,4) + 2*C_fs/(gamma-1);
Rnim = u(round(0.25*I_MAX):round(0.75*I_MAX),J_MAX-1,3) + ...
     +u(round(0.25*I_MAX):round(0.75*I_MAX), J_MAX-1,4) - 2*C_fs/(gamma-1);
%unbi = (Rnbp + Rnim)/2;
cbi = (Rnbp - Rnim)*0.25*(gamma-1);

P = (gamma -1)*(U(:,:,4) - 0.5*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
uvec = ones(I_MAX,J_MAX,2);
uvec(:,:,1) = u_fs; 
uvec(:,:,2) = v_fs;
nu = zeros(I_MAX,1);
E2 = zeros(I_MAX-1,J_MAX-1,4);
E4 = zeros(I_MAX-1,J_MAX-1,4);
k2 = 1/4; 
k4 = 1/256;

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

for nind = 0:3
    for ii = 1:I_MAX-1
        for jj = 1:J_MAX-1
            magds = sqrt(cell_n(ii,jj,(1+nind*2))^2 + cell_n(ii,jj,(2+nind*2))^2);
            if ii == 1
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + C_fs*magds)*max([nu(I_MAX),nu(ii), nu(ii+1), nu(ii+2)]);
            elseif ii == I_MAX-1
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + C_fs*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(1)]);
            else
                E2(ii,jj,nind+1) = 0.5*k2*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                    + C_fs*magds)*max([nu(ii-1),nu(ii), nu(ii+1), nu(ii+2)]);
            end

            E4(ii,jj,nind+1) = max([0,(0.5*k4*(dot(uvec(ii,jj,:),cell_n(ii,jj,(1+nind*2):(2+nind*2)))...
                + C_fs*magds) - E2(ii,jj,nind+1))]);

            if ii== 1
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(ii+2,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:)...
                    + 3*U(ii,jj,nind+1,:) - U(I_MAX-1,jj,nind+1,:)));
            elseif ii == I_MAX-2
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(1,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:) + ...
                    3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));

            elseif ii == I_MAX-1
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(1,jj,nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(2,jj,nind+1,:) - 3*U(1,jj,nind+1,:) + ...
                    3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));
            else
                D(ii,jj,nind+1,:) = E2(ii,jj,nind+1).*(U(ii+1,jj, nind+1,:) - U(ii,jj,nind+1,:)) - ...
                    (E4(ii,jj,nind+1).*(U(ii+2,jj,nind+1,:) - 3*U(ii+1,jj,nind+1,:)...
                    + 3*U(ii,jj,nind+1,:) - U(ii-1,jj,nind+1,:)));
            end
        end
    end
end


% %from here assuming the normal boundary conditions
% U = U_IC; G = G_IC; F = F_IC;
% %calculate Reimann invariants based on freestream
% Rnbpi = u_IC(0.25*I_MAX:0.75*I_MAX,J_MAX-1,3) + ...
%     +u_IC(0.25*I_MAX:0.75*I_MAX, J_MAX-1,4) + 2*cfs/(gamma-1);
% Rnim = u_IC(0.25*I_MAX:0.75*I_MAX,J_MAX-1,3) + ...
%     +u_IC(0.25*I_MAX:0.75*I_MAX, J_MAX-1,4) - 2*cfs/(gamma-1);
% %unbi = (Rnbp + Rnim)/2;
% cbi = (Rnbp - Rnim)*0.25*(gamma-1);
% %use streamline constants to calculate interior values
% Hfs = (gamma/(gamma-1))*Pfs/rhofs + 0.5*ufs^2;


end

