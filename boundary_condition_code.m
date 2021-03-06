function UB2 = boundary_condition_code(I_MAX, J_MAX, U,S)



%setting up the boundary conditions 
rho_fs = 0.4135; u_fs = 258.4; v_fs = 0.0; C_fs = 304.025; M_fs = 0.85;
p_fs = 27300.0;
E_fs = (1/(1.4-1))*(p_fs/rho_fs) + (u_fs^2 + v_fs^2)/2;
gamma = 1.4; 
P = (gamma -1)*(U(:,:,4) - 0.5*(U(:,:,2).^2 + U(:,:,3).^2)./U(:,:,1));
C = sqrt(gamma*P./U(:,:,1));
u_i = U((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),J_MAX-1,2)./U((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,1);
ibnorm = sqrt(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1).^2 + ...
    S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2).^2);
unfs = -u_fs*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm);
Rnbpi = unfs + ones((0.5*(I_MAX-1)),1)*2*C_fs/(gamma-1);

u_in =  -u_i.*((S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)+...
    S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2))./ibnorm);

Rnimi = u_in - 2*C((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),J_MAX-1)/(gamma-1);

unbi = (Rnbpi + Rnimi)/2;
cbi = ((gamma-1)/4)*(Rnbpi + Rnimi);

%tangential velocities
ulbi = u_fs*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm);
%pressure vector at boundary set to interior adjacent cell
Pbi = P((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)), J_MAX-1);
%solved for rho based on that
rhobi = gamma*Pbi./(cbi.^2);
%set the x and y velocities for each of the components
ubi = ulbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm) - ...
    unbi.*(S(0.25*(I_MAX-1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm);
vbi = -ulbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,1)./ibnorm) - ...
    unbi.*(S((0.25*(I_MAX-1)+1):0.75*(I_MAX-1),J_MAX-1,4,2)./ibnorm);
Ebi = (1/(gamma-1))*(Pbi./rhobi) + 0.5*(ubi.^2 + vbi.^2);
%put together U
Ubi = zeros((0.5*I_MAX-1),4);


Ubi(:,1) = rhobi;
Ubi(:,2) = rhobi.*ubi;
Ubi(:,3) = rhobi.*vbi;
Ubi(:,4) = rhobi.*Ebi;

% outer boundary conditions, upper quadrant
u_ou = U(((0.75*(I_MAX-1)+1)):(I_MAX-1),J_MAX-1,2)./U((0.75*(I_MAX-1)):(I_MAX-1),J_MAX-1,1);
obnormu = sqrt(S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1).^2 ...
+ S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2).^2);
unfsou = u_fs*(S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1)./obnormu);
Rnbpou = -abs(unfsou) + ones((0.25*(I_MAX-1)),1)*P_fs/(rho_fs*C_fs);
%normal velocity
u_onu =  u_ou.*((S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1)+...
    S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2))./obnormu);
%tangential velocity
u_olu =  u_ou.*((-S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2)+...
    S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1))./obnormu);
Rnimou = u_onu - 2*C((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1)/(gamma-1);

unbou = -(Rnbpou + Rnimou)/2; %note this is an absolute value/magnitude
cbou = -abs(unbou) - Rnbpou;

%outer boundary conditions, lower quadrant 
u_ol = U(1:(0.25*(I_MAX-1)),J_MAX-1,2)./U(1:(0.25*(I_MAX-1)),J_MAX-1,1);
obnorml = sqrt(S(1:(0.25*(I_MAX-1)),J_MAX-1,7).^2 ...
+ S(1:(0.25*(I_MAX-1)),J_MAX-1,8).^2);
unfsol = ufs*(S(1:(0.25*(I_MAX-1)),J_MAX-1,7)./obnorml);
Rnbpol = -abs(unfsol) + ones((0.25*(I_MAX-1)),1)*P_fs/(rho_fs*C_fs);

u_onl =  u_ol.*((S(1:(0.25*(I_MAX-1)),J_MAX-1,7)+...
    S(1:(0.25*(I_MAX-1)),J_MAX-1,8))./obnorml);
%tangential velocity lower boundary
u_oll =  u_ol.*((-S(1:(0.25*(I_MAX-1)),J_MAX-1,8)+...
    S(1:(0.25*(I_MAX-1)),J_MAX-1,7))./obnorml);

Rnimol = u_onl - 2*C(1:(0.25*(I_MAX-1)),J_MAX-1)/(gamma-1);

unbol = -(Rnbpol + Rnimol)/2; %note this is an absolute value/magnitude
cbol = -abs(unbol) - Rnbpol;


unbo = [unbou unbol];
ulbo = [u_olu u_oll];
cbo = [cbou cbol];
obnorm = [obnormo obnorml];
%make the vectors of the normal/tangential vector components
Snormx = [S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,1) S(1:(0.25*(I_MAX-1)),J_MAX-1,4,1)];
Snormy = [S((0.75*(I_MAX-1)+1):(I_MAX-1),J_MAX-1,4,2) S(1:(0.25*(I_MAX-1)),J_MAX-1,4,2)];

%pressure vector at boundary set to interior adjacent cell
Pbo = [P((0.75*(I_MAX-1)+1):I_MAX-1, J_MAX-1) P(1:(0.25*(I_MAX-1)), J_MAX-1)];
%solved for rho based on that
rhobo = gamma*Pbo./(cbo.^2);
%set the x and y velocities for each of the components
ubo = -ulbo.*(Snormy./obnorm) + unbo.*(Snormx./obnorm);
vbo = ulbo.*(Snormx./obnorm) + unbo.*(Snormy./obnorm);
Ebo = (1/(gamma-1))*(Pbo./rhobo) + 0.5*(ubo.^2 + vbo.^2);
%put together U, F, G
Ubo = zeros((0.5*I_MAX-1),4);

Ubo(:,1) = rhobo;
Ubo(:,2) = rhobo.*ubo;
Ubo(:,3) = rhobo.*vbo;
Ubo(:,4) = rhobo.*Ebo;

UB2 = zeros(I_MAX-1,4);
UB2(1:(0.25*(I_MAX-1)),:) = Ubo((0.25*(I_MAX-1)+1):(0.5*(I_MAX-1)),:);
UB2((0.25*(I_MAX-1)+1):(0.75*(I_MAX-1)),:) = Ubi;
UB2((0.75*(I_MAX-1)+1):I_MAX-1,:) = Ubo(1:(0.25*(I_MAX-1)),:);

end