function UB2 = Compute_UB2(I_MAX, J_MAX, U,S)



%setting up the boundary conditions 
rho_fs = 0.4135; u_fs = 258.4; v_fs = 0.0; C_fs = 304.025; M_fs = 0.85;
p_fs = 27300.0;
gamma = 1.4; 
E_fs = (1/(1.4-1))*(p_fs/rho_fs) + (u_fs^2 + v_fs^2)/2;
Hfs = E_fs + p_fs/rho_fs;
Sfs = p_fs/(rho_fs^gamma);


P = (gamma -1)*(U(:,:,4) - 0.5*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
C = sqrt(gamma*P./U(:,:,1));
u_i = U((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),J_MAX-1,2)./U((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,1);
v_i = U((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),J_MAX-1,3)./U((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,1);
ibnorm = sqrt(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1).^2 + ...
    S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2).^2);
unfs = -u_fs*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm);
Rnbpi = unfs + ones((0.5*(I_MAX-1)),1)*2*C_fs/(gamma-1);

u_in =  (u_i.*(-S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)))./ibnorm+...
    (v_i.*(-S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)))./ibnorm;

Rnimi = u_in - 2*C((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),J_MAX-1)/(gamma-1);

unbi = (Rnbpi + Rnimi)/2;
cbi = ((gamma-1)/4)*(Rnbpi - Rnimi);

%tangential velocities
ulbi = u_fs*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm);
%pressure vector at boundary set to interior adjacent cell
%rhobi = U((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)), J_MAX-1,1);
%solved for rho based on entropy
rhobi = ((cbi.^2)/(gamma*Sfs)).^(1/(gamma-1));
%set the x and y velocities for each of the components
ubi = ulbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm) - ...
    unbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm);
vbi = -ulbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm) - ...
    unbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm);

%Pbi = (gamma/(gamma-1))*(Hfs - 0.5*(ubi.^2 + vbi.^2)).*rhobi;
Pbi = rhobi.*(cbi.^2)/gamma;
Ebi = (1/(gamma-1))*(Pbi./rhobi) + 0.5*(ubi.^2 + vbi.^2);
%put together U
Ubi = zeros((0.5*(I_MAX-1)),4);


Ubi(:,1) = rhobi;
Ubi(:,2) = rhobi.*ubi;
Ubi(:,3) = rhobi.*vbi;
Ubi(:,4) = rhobi.*Ebi;

% outer boundary conditions, upper quadrant
u_ou = U(((0.75*(I_MAX-1)+1)):(I_MAX-1),J_MAX-1,2)./U((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,1);
v_ou = U(((0.75*(I_MAX-1)+1)):(I_MAX-1),J_MAX-1,3)./U((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,1);
obnormu = sqrt(S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1).^2 ...
+ S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2).^2);
unfsou = u_fs*(S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1)./obnormu);
Rnbpou = -abs(unfsou) - 2*C_fs/(gamma-1);
%normal velocity
u_onu =  (u_ou.*S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1))./obnormu+...
    (v_ou.*S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2))./obnormu;
%tangential velocity
u_olu =  (u_ou.*-S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2))./obnormu+...
    (v_ou.*S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1))./obnormu;

Rnimou = -abs(u_onu) + 2*C((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1)/(gamma-1);

unbou = -(Rnbpou + Rnimou)/2; %note this is an absolute value/magnitude
cbou = ((gamma-1)/4)*(Rnbpou - Rnimou);

%outer boundary conditions, lower quadrant 
u_ol = U(1:(0.25*(I_MAX-1)),J_MAX-1,2)./U(1:(0.25*(I_MAX-1)),J_MAX-1,1);
v_ol = U(1:(0.25*(I_MAX-1)),J_MAX-1,3)./U(1:(0.25*(I_MAX-1)),J_MAX-1,1);
obnorml = sqrt(S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1).^2 ...
+ S(1:(0.25*(I_MAX-1)),J_MAX-1,4,2).^2);
unfsol = u_fs*(S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1)./obnorml);
Rnbpol = -abs(unfsol)  - 2*C_fs/(gamma-1);

u_onl =  (u_ol.*S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1))./obnorml+...
    (v_ol.*S(1:(0.25*(I_MAX-1)),J_MAX-1,4,2))./obnorml;
%tangential velocity
u_oll =  (u_ol.*-S(1:(0.25*(I_MAX-1)),J_MAX-1,4,2))./obnorml+...
    (v_ol.*S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1))./obnorml;

Rnimol = -abs(u_onl) + 2*C(1:(0.25*(I_MAX-1)),J_MAX-1)/(gamma-1);

unbol = -(Rnbpol + Rnimol)/2; %note this is an absolute value/magnitude
cbol = ((gamma-1)/4)*(Rnbpol - Rnimol);


unbo = [unbou' unbol'];
unbo = unbo';
ulbo = [u_olu' u_oll'];
ulbo = ulbo';
cbo = [cbou' cbol'];
cbo = cbo';
obnorm = [obnormu' obnorml'];
obnorm = obnorm';
%make the vectors of the normal/tangential vector components
Snormx = [S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1)' S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1)'];
Snormy = [S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2)' S(1:(0.25*(I_MAX-1)),J_MAX-1,4,2)'];
Snormx = Snormx';
Snormy = Snormy';

%pressure vector at boundary set to interior adjacent cell
%rhobo = [U((0.75*(I_MAX-1)+1):I_MAX-1, J_MAX-1,1)' U(1:(0.25*(I_MAX-1)), J_MAX-1,1)'];
%rhobo = rhobo';
rhobo = ((cbo.^2)/(gamma*Sfs)).^(1/(gamma-1));
%solved for rho based on that

%set the x and y velocities for each of the components
ubo = -ulbo.*(Snormy./obnorm) + unbo.*(Snormx./obnorm);
vbo = ulbo.*(Snormx./obnorm) + unbo.*(Snormy./obnorm);
%Pbo = (gamma/(gamma-1))*(Hfs - 0.5*(ubo.^2 + vbo.^2)).*rhobo;
Pbo = rhobo.*(cbo.^2)/gamma;
Ebo = (1/(gamma-1))*(Pbo./rhobo) + 0.5*(ubo.^2 + vbo.^2);
%put together U, F, G
Ubo = zeros((0.5*(I_MAX-1)),4);

Ubo(:,1) = rhobo;
Ubo(:,2) = rhobo.*ubo;
Ubo(:,3) = rhobo.*vbo;
Ubo(:,4) = rhobo.*Ebo;

UB2 = zeros(I_MAX-1,4);
UB2(1:(0.25*(I_MAX-1)),:) = Ubo((0.25*(I_MAX-1)+1):(0.5*(I_MAX-1)),:);
UB2((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),:) = Ubi;
UB2((0.75*(I_MAX-1)+1):I_MAX-1,:) = Ubo(1:(0.25*(I_MAX-1)),:);

end