function UB1 = Compute_UB1( I_MAX, J_MAX, U,S,X)

%extrapolate each quantity
gamma = 1.4;
P = (gamma -1)*(U(:,:,4) - 0.5*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
C = sqrt(gamma*P./U(:,:,1));
E = U(:,:,4)./U(:,:,1);
dist = sqrt((X(:,2,2) - X(:,1,2)).^2 + (X(:,2,1) - X(:,1,1)).^2);
P_slopes = (P(:,2) - P(:,1))./dist;
rhoslopes = (U(:,2,1) - U(:,1,1))./dist;
Eslopes = (E(:,2) - E(:,1))./dist;

u = U(:,1,2)./U(:,1,1);
v = U(:,1,3)./U(:,1,1);
UB1 = zeros(I_MAX-1,4);
boundnorms = sqrt(S(:,1,3,2).^2 + S(:,1,3,1).^2);
%tangential quantities for the surface
ul = u.*S(:,1,3,2)./boundnorms - v.*S(:,1,3,1)./boundnorms; 
% now recompute the u and v for physical 
uba = ul.*S(:,1,3,2)./boundnorms;
vba = -ul.*S(:,1,3,1)./boundnorms;

%ratio of the distances to find the next distance for extrapolation
dist2 = sqrt((X(:,3,2) - X(:,2,2)).^2 + (X(:,3,1) - X(:,2,1)).^2);
distrat = dist./dist2;
%theta = arctan(X(:,2,2)./X(:,2,1));
distin = distrat.*dist;
%xin = distin.*cos(theta);
%yin = distin.*sin(theta);

%do the extrapolation of values 
Pba = P(:,1) + ((distin - dist)./(dist2-dist)) .* (P(:,2) - P(:,1));
rhoba = U(:,1,1) + ((distin - dist)./(dist2-dist)) .* (U(:,2,1) - U(:,1,1));
Eba  = E(:,1) + ((distin - dist)./(dist2-dist)) .* (E(:,2) - E(:,1));

UB1(:,1) = rhoba;
UB1(:,2) = rhoba.*uba;
UB1(:,3) = rhoba.*vba;
UB1(:,4) = rhoba.*Eba;

end