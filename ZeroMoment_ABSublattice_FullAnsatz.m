close all
clearvars

global  Nu Nk 
global  Mekx Meky Mweightk
global  MJx MJy MJz T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Brillouin zone cutting setting

N1 = 36; % number of points along k_1-axis (oblique coordinate) % 144 % 7
N2 = N1;
Nu = N1*N2; % number of unit cell
Ns = 2*Nu; % number of Lattice sites


if mod(N1,2) == 0

[ky_grid, kx_grid] = meshgrid(2:round(N2/2), 1:N1);
kx_BZ = 2 * (kx_grid(:) - 1) / N1 * pi;
ky_BZ = 2 * (ky_grid(:) - 1) / N2 * pi;


kk_edge = (2:round(N2/2));
kx_edge = [2*(kk_edge-1)/N1*pi, 2*(kk_edge-1)/N1*pi];
ky_edge = [zeros(1,length(kk_edge)), pi*ones(1,length(kk_edge))];

kx0 = [kx_BZ.', kx_edge];
ky0 = [ky_BZ.', ky_edge];
NBz = length(kx0);

% Unified all k-points: first NBz have weight 2,  
% 4 particle-hole symmetric points have weight 1 :
% -k = k + 2*n*pi => k=(0,0),(0,pi),(pi,0),(pi,pi);
Mekx = [exp(-1i*kx0),1, 1, -1, -1];
Meky = [exp(-1i*ky0),1, -1, 1, -1];
Mweightk = [2*ones(1,NBz), ones(1,4)];
Nk = NBz+4;

else

[ky_grid, kx_grid] = meshgrid(2:round( (N2+1)/2 ), 1:N1);
kx_BZ = 2 * (kx_grid(:) - 1) / N1 * pi;
ky_BZ = 2 * (ky_grid(:) - 1) / N2 * pi;

kk_edge = (2:round( (N2+1)/2 ));
kx_edge = 2*(kk_edge-1)/N1*pi;
ky_edge = zeros(1,length(kk_edge));

kx0 = [kx_BZ.', kx_edge];
ky0 = [ky_BZ.', ky_edge];
NBz = length(kx0);

% Unified all k-points: first NBz have weight 2,  
% 1 particle-hole symmetric points have weight 1. k=(0,0);
Mekx = [exp(-1i*kx0),1];
Meky = [exp(-1i*ky0),1];
Mweightk = [2*ones(1,NBz), 1];
Nk = NBz+1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solved QSL MF Ansatz 

rr0 = [0.3	-0.524863394	0	0.080859055	0	-0.662396223	0	1]; % 2 cones
% rr0 = [0.3	-0.509863986	-0.098167284	0.754200435	0.084402609	-0.276428067	-0.144750182	1]; % 8 cones
% rr0 = [0.3	-0.521998077	0.054532335	0.825567732	-0.033218347	-0.262831995	0.04692293	1 ]; % 14_1 cones
% rr0 = [0.3	-0.341432501	-0.004289107	0.370827978	-0.260717961	-0.35438937	0.258113385	1]; % 14_2 cones
% rr0 = [0.3	-0.419865341	0.007071734	0.405669368	0.205573093	-0.413253422	-0.218953425	1]; % 8+8 cones
% rr0 = [0.3	-0.285097553	-0.00968513	0.304467813	-0.299625248	-0.294782683	0.289940118	1]; % 14+2 cones
% rr0 = [0.3	-0.524863394	0	0.847804625	0	-0.248934499	0	1]; % 20_1 cones
% rr0 = [0.3	-0.407860725	0.053371595	0.252725455	0.18737562	-0.542058344	-0.150538265	1]; % 20_2 cones
% rr0 = [0.3	-0.524863396	0	0.572639563	0	-0.459278014	0	1]; % 20_3 cones
% rr0 = [0.3	-0.440881357	-0.027323882	0.350925451	-0.14964179	-0.505132307	0.162800282	1]; % 20_4 cones
% rr0 = [0.3	-0.459729283	0.046615123	0.568111446	-0.155095844	-0.380883429	0.166746217	1]; % 26 cones
% rr0 = [0.3	0.017147549	-0.017147549	0.017147549	0.357284377	0	-0.357284377	1]; % 32 cones
% rr0 = [-0.6	-0.519463323	-0.042690589	0.61284965	-0.068167345	0.472734569	0.021380137	1];

% rr0 = [0.0  -0.528240628	-4.01E-10	0	9.82E-10	0	-6.10E-10	1]; %pure K = -1;
% rr0 = [0.0  0.528240628	-4.01E-10	0	9.82E-10	0	-6.10E-10	1]; %pure K = 1;

% rr0=[ 0.0, 2*rand(1,6)-1 ]; rand trial

chiz_cc = rr0(2);
chiz_xx = rr0(3);
chiz_zz = rr0(4);
chiz_cx = rr0(5);
chiz_xy = rr0(6);
chiz_xz = rr0(7);

chi0_z = [chiz_cc, chiz_cx, -chiz_cx,       0;
          chiz_cx, chiz_xx, chiz_xy , chiz_xz;
         -chiz_cx, chiz_xy, chiz_xx , chiz_xz;
                0, chiz_xz, chiz_xz , chiz_zz]; 

Mcyc = [1,0,0,0;[0;0;0],[0 1 0;0 0 1;1 0 0]];
chi0_x = Mcyc'*chi0_z*Mcyc;
chi0_y = Mcyc*chi0_z*Mcyc';

MF_sol0 = [chi0_x(:).', chi0_y(:).', chi0_z(:).', zeros(1,6)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('TolFun',1e-8);
has_previous_solution = false;

ww = 0.3;
% www = (0.3:0.01:0.6);

split         = length(ww);
Data0_full    = NaN(split,54+2);
Data0_EG      = NaN(split,6+1);
Data0_DD      = NaN(split,3+1);
Data0_MFCheck = NaN(split,5+1);

sii=0;

for ii = 1:split

J = 0;
K = -1;
%Ga = 0.1;
Ga_p = 0;
T = 1e-6;

%J = ww(ii);
%K = ww(ii);
Ga = ww(ii);
%Ga_p = ww(ii);
%T = ww(ii); 

% RJ=[1,  4/9, -4/9,  4/9;
%     0, -1/3,  4/3, -4/3;
%     0   4/9,  5/9,  4/9;
%      0 -2/9,  2/9,  7/9];
% LJ=[J;K;Ga;Ga_p];
% LLJ=RJ*LJ;
% 
% J = LLJ(1);
% K = LLJ(2);
% Ga = LLJ(3);
% Ga_p = LLJ(4); hidden paramters transformation 

MJx = [ J+K, Ga_p, Ga_p;
       Ga_p,    J,   Ga;
       Ga_p,   Ga,    J]/4; % 1/4 : s*s = 1/4*sigma*sigma
MJy = [   J, Ga_p,   Ga;
       Ga_p,  J+K, Ga_p;
         Ga, Ga_p,    J]/4;
MJz = [   J,   Ga, Ga_p;
         Ga,    J, Ga_p;
       Ga_p, Ga_p,  J+K]/4;

if has_previous_solution == false
    v_init = MF_sol0;
else
    v_init = v_previous;
end
   
[v,fval,exitflag] = fsolve(@ MF_Eq, v_init, options);

if exitflag == 1
    has_previous_solution = true;
    v_previous = v;
    sii=sii+1;
end 
 
if exitflag == 0
    break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data0_full(ii,:) = [ww(ii), v, exitflag];

chi_x = reshape(v(1:16),[4,4]); 
chi_y = reshape(v(17:32),[4,4]);
chi_z = reshape(v(33:48),[4,4]); 

DD_x = det(chi_x); DD_y = det(chi_y); DD_z = det(chi_z);
Data0_DD(ii,1) = ww(ii);
Data0_DD(ii,2) = DD_x;
Data0_DD(ii,3) = DD_y;
Data0_DD(ii,4) = DD_z;

obs = obs_Eq(v);       
Ev0 = (1/2)*(obs.gJep_x+obs.gJep_y+obs.gJep_z);
EvT = Ev0+(obs.F0-obs.E0)/Ns; % RMFT energy at finite T
E1  = 0.5*sum(obs.tx.*chi_x+obs.ty.*chi_y+obs.tz.*chi_z,'all');
Data0_EG(ii,:) = [ww(ii), EvT, Ev0, obs.S0/Ns,...
                 exitflag, sqrt(mean(fval.^2)), obs.E0/Ns-E1]; % if MF converged then E0 = E1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data0_MFCheck(ii,1) = ww(ii);
Data0_MFCheck(ii,2) = norm(obs.Sm_A)+norm(obs.Sm_B);
Data0_MFCheck(ii,3) = norm(obs.Sq_A)+norm(obs.Sq_B);
Data0_MFCheck(ii,4) = norm(chi_x - obs.IS_x,'fro');
Data0_MFCheck(ii,5) = norm(chi_y - obs.IS_y,'fro');
Data0_MFCheck(ii,6) = norm(chi_z - obs.IS_z,'fro');

end

Data_full    = Data0_full(1:sii,:);
Data_EG      = Data0_EG(1:sii,:);
Data_DD      = Data0_DD(1:sii,:);
Data_MFCheck = Data0_MFCheck(1:sii,:);

% 
if ww(split)<ww(1)
Data_full    = fliplr(Data_full.').';
Data_EG      = fliplr(Data_EG.').';
Data_DD      = fliplr(Data_DD.').';
Data_MFCheck = fliplr(Data_MFCheck.').';
end


function F = MF_Eq( v )

global  Nu Nk 
global  Mekx Meky Mweightk
global  MJx MJy MJz T

% chi_x = [v(1), v(5),  v(9), v(13),
%         v(2), v(6), v(10), v(14),
%         v(3), v(7), v(11), v(15),
%         v(4), v(8), v(12), v(16)]  

chi_x = reshape(v(1:16),[4,4]); 
chi_y = reshape(v(17:32),[4,4]);
chi_z = reshape(v(33:48),[4,4]);
La_A = v(49:51); 
La_B = v(52:54);

[Cx,tx0] = build_block(chi_x,MJx); 
[Cy,ty0] = build_block(chi_y,MJy);
[Cz,tz0] = build_block(chi_z,MJz);

ep_x = sum(MJx.*Cx,'all');
ep_y = sum(MJy.*Cy,'all');
ep_z = sum(MJz.*Cz,'all');

DD_x = det(chi_x); % <D_i D_j>_x
DD_y = det(chi_y);
DD_z = det(chi_z);          
DD_s = DD_x+DD_y+DD_z;

gJ_x = 4/(1+2*DD_s-DD_x);
gJ_y = 4/(1+2*DD_s-DD_y);
gJ_z = 4/(1+2*DD_s-DD_z);

dDD_x = CoMatrix4(chi_x); % \partial <D_i D_j>_x  /  \partial \chi^{\mu\nu}
dDD_y = CoMatrix4(chi_y);
dDD_z = CoMatrix4(chi_z);

gJep_x=gJ_x*ep_x;
gJep_y=gJ_y*ep_y;
gJep_z=gJ_z*ep_z;
      
tx = gJ_x*tx0-(gJ_x*gJep_x+2*gJ_y*gJep_y+2*gJ_z*gJep_z)/4*dDD_x;
ty = gJ_y*ty0-(gJ_y*gJep_y+2*gJ_z*gJep_z+2*gJ_x*gJep_x)/4*dDD_y;
tz = gJ_z*tz0-(gJ_z*gJep_z+2*gJ_x*gJep_x+2*gJ_y*gJep_y)/4*dDD_z; 
      
Htx = 0.5*1i*tx; % 0.5: ib_A*b_B = (ib_A*b_B-ib_B*b_A)/2 => ib_A*b_B -> 0.5*(Htx + Htx')
Hty = 0.5*1i*ty;
Htz = 0.5*1i*tz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heff_A \cdot \sigma_A +  lambda_A \cdot Q_A
% heff_A = \partial E_v  /  \partial \vec{m}_A
% = -h_ext/2 + (gJx*MJx+...)*mB +...

heff_A = zeros(1,3); 
heff_B = zeros(1,3);

HAA = (1i*[0, -heff_A(1), -heff_A(2), -heff_A(3);
           0,          0, -heff_A(3),          0;
           0,          0,          0, -heff_A(1);
           0, -heff_A(2),          0,          0]/4)...
     +(1i*[0,    La_A(1),    La_A(2),    La_A(3);
           0,          0,   -La_A(3),          0;
           0,          0,          0,   -La_A(1);
           0,   -La_A(2),          0,          0]/4);
HAA = HAA+HAA';


HBB = (1i*[0, -heff_B(1), -heff_B(2), -heff_B(3);
           0,          0, -heff_B(3),          0;
           0,          0,          0, -heff_B(1);
           0, -heff_B(2),          0,          0]/4)...
     +(1i*[0,    La_B(1),    La_B(2),    La_B(3);
           0,          0,   -La_B(3),          0;
           0,          0,          0,   -La_B(1);
           0,   -La_B(2),          0,          0]/4);
HBB = HBB+HBB';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Schi_x = zeros(4);
Schi_y = zeros(4);
Schi_z = zeros(4);
SG_A = zeros(4);
SG_B = zeros(4);
        
  for ii = 1:Nk
       
        ekx = Mekx(ii);
        eky = Meky(ii);
        wk  = Mweightk(ii);

        HAB = Htx*ekx+Hty*eky+Htz;        
        H = 4*[HAA,HAB;HAB',HBB]; % 4:  E_k = 4*lambda_k
        
        [VV,DD] = eig(H);          
        EE = diag(DD).';       
        fermi = 0.5*(1-tanh(EE/(2*T)));
            
        WW = conj(VV);
        AvgAB = wk*(WW(1:4, :) .* fermi) * (VV(5:8, :).'); % <b^{\mu}_A *b^{\nu}_B}>

        Schi_x = Schi_x+ekx*AvgAB;
        Schi_y = Schi_y+eky*AvgAB;
        Schi_z = Schi_z+AvgAB;              
  
        SG_A = SG_A+wk*(WW(1:4, :) .* fermi) * (VV(1:4, :).');
        SG_B = SG_B+wk*(WW(5:8, :) .* fermi) * (VV(5:8, :).');        
    
  end

IS_x = -2/Nu*imag(Schi_x);
IS_y = -2/Nu*imag(Schi_y);
IS_z = -2/Nu*imag(Schi_z);

IG_A = -2/Nu*imag(SG_A);
IG_B = -2/Nu*imag(SG_B);

WSchi_x = (IS_x(:)).';
WSchi_y = (IS_y(:)).';
WSchi_z = (IS_z(:)).';
  
% Sm_A = ([IGA(2,1)-IGA(3,4), ...
%          IGA(3,1)-IGA(4,2), ...
%          IGA(4,1)-IGA(2,3)])/2;
Sq_A = -([IG_A(2,1)+IG_A(3,4), ...
          IG_A(3,1)+IG_A(4,2), ...
          IG_A(4,1)+IG_A(2,3)])/2;
% Sm_B = ([IGB(2,1)-IGB(3,4), ...
%          IGB(3,1)-IGB(4,2), ...
%          IGB(4,1)-IGB(2,3)])/2;
Sq_B = -([IG_B(2,1)+IG_B(3,4), ...
          IG_B(3,1)+IG_B(4,2), ...
          IG_B(4,1)+IG_B(2,3)])/2;

F = zeros(1,54);
F(1:16)  = v(1:16)-WSchi_x;
F(17:32) = v(17:32)-WSchi_y;
F(33:48) = v(33:48)-WSchi_z;
F(49:51) = Sq_A;
F(52:54) = Sq_B;

end

function obs = obs_Eq( v )

global  Nu Nk 
global  Mekx Meky Mweightk
global  MJx MJy MJz T

chi_x = reshape(v(1:16),[4,4]);
chi_y = reshape(v(17:32),[4,4]);
chi_z = reshape(v(33:48),[4,4]);
La_A = v(49:51); 
La_B = v(52:54); 

[Cx,tx0] = build_block(chi_x,MJx); 
[Cy,ty0] = build_block(chi_y,MJy);
[Cz,tz0] = build_block(chi_z,MJz);

ep_x = sum(MJx.*Cx,'all');
ep_y = sum(MJy.*Cy,'all');
ep_z = sum(MJz.*Cz,'all');

DD_x = det(chi_x); % <D_i D_j>_x
DD_y = det(chi_y);
DD_z = det(chi_z);          
DD_s = DD_x+DD_y+DD_z;
gJ_x = 4/(1+2*DD_s-DD_x);
gJ_y = 4/(1+2*DD_s-DD_y);
gJ_z = 4/(1+2*DD_s-DD_z);

dDD_x = CoMatrix4(chi_x); % \partial <D_i D_j>_x  /  \partial \chi^{\mu\nu}
dDD_y = CoMatrix4(chi_y);
dDD_z = CoMatrix4(chi_z);

gJep_x = gJ_x*ep_x; obs.gJep_x = gJep_x;
gJep_y = gJ_y*ep_y; obs.gJep_y = gJep_y;
gJep_z = gJ_z*ep_z; obs.gJep_z = gJep_z;

      
tx = gJ_x*tx0-(gJ_x*gJep_x+2*gJ_y*gJep_y+2*gJ_z*gJep_z)/4*dDD_x; obs.tx = tx;
ty = gJ_y*ty0-(gJ_y*gJep_y+2*gJ_z*gJep_z+2*gJ_x*gJep_x)/4*dDD_y; obs.ty = ty;
tz = gJ_z*tz0-(gJ_z*gJep_z+2*gJ_x*gJep_x+2*gJ_y*gJep_y)/4*dDD_z; obs.tz = tz;
      
Htx = 0.5*1i*tx; % 0.5: ib_A*b_B = (ib_A*b_B-ib_B*b_A)/2 => ib_A*b_B -> 0.5*(Htx + Htx')
Hty = 0.5*1i*ty;
Htz = 0.5*1i*tz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heff_A \cdot \sigma_A +  lambda_A \cdot Q_A
% heff_A = \partial E_v  /  \partial \vec{m}_A
% = -h_ext/2 + (gJx*MJx+...)*mB +...

heff_A = zeros(1,3); 
heff_B = zeros(1,3);

HAA = (1i*[0, -heff_A(1), -heff_A(2), -heff_A(3);
           0,          0, -heff_A(3),          0;
           0,          0,          0, -heff_A(1);
           0, -heff_A(2),          0,          0]/4)...
     +(1i*[0,    La_A(1),    La_A(2),    La_A(3);
           0,          0,   -La_A(3),          0;
           0,          0,          0,   -La_A(1);
           0,   -La_A(2),          0,          0]/4);
HAA = HAA+HAA';


HBB = (1i*[0, -heff_B(1), -heff_B(2), -heff_B(3);
           0,          0, -heff_B(3),          0;
           0,          0,          0, -heff_B(1);
           0, -heff_B(2),          0,          0]/4)...
     +(1i*[0,    La_B(1),    La_B(2),    La_B(3);
           0,          0,   -La_B(3),          0;
           0,          0,          0,   -La_B(1);
           0,   -La_B(2),          0,          0]/4);
HBB = HBB+HBB';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Schi_x = zeros(4);
Schi_y = zeros(4);
Schi_z = zeros(4);
SG_A = zeros(4);
SG_B = zeros(4);

F0 = 0;
S0 = 0;
E0 = 0;
        
  for ii = 1:Nk
       
        ekx = Mekx(ii);
        eky = Meky(ii);
        wk  = Mweightk(ii);

        HAB = Htx*ekx+Hty*eky+Htz;
        
        H = 4*[HAA,HAB;HAB',HBB]; % 4:  E_k = 4*lambda_k
        
        [VV,DD] = eig(H);          
        EE = diag(DD).';
        
        fermi = 0.5*(1-tanh(EE/(2*T)));
        
        AEE = abs(EE);        
        fermiA = 0.5*(1-tanh(AEE/(2*T)));  
        EE0 = wk/2*sum(AEE)/2;
        F0 = F0-wk/2*T*sum(log(1+exp(-AEE/T)))-EE0;
        E0 = E0+wk/2*sum(AEE.*fermiA)-EE0;
%       S0 = S0-sum((1-fermiA).*log(1-fermiA)+fermiA.*log(fermiA));
        S0 = S0+wk/2*sum(AEE/T./(1+exp(AEE/T))+log(1+exp(-AEE/T)));
          
        WW = conj(VV);
        AvgAB = (WW(1:4, :) .* fermi) * (VV(5:8, :).'); % <b^{\mu}_A *b^{\nu}_B}>

        Schi_x = Schi_x+wk*ekx*AvgAB;
        Schi_y = Schi_y+wk*eky*AvgAB;
        Schi_z = Schi_z+wk*AvgAB;              
  
        SG_A = SG_A+wk*(WW(1:4, :) .* fermi) * (VV(1:4, :).');
        SG_B = SG_B+wk*(WW(5:8, :) .* fermi) * (VV(5:8, :).');        
    
  end

obs.IS_x = -2/Nu*imag(Schi_x);
obs.IS_y = -2/Nu*imag(Schi_y);
obs.IS_z = -2/Nu*imag(Schi_z);

IG_A = -2/Nu*imag(SG_A);
IG_B = -2/Nu*imag(SG_B);
  
obs.Sm_A = ([IG_A(2,1)-IG_A(3,4), ...
             IG_A(3,1)-IG_A(4,2), ...
             IG_A(4,1)-IG_A(2,3)])/2;
obs.Sq_A = -([IG_A(2,1)+IG_A(3,4), ...
              IG_A(3,1)+IG_A(4,2), ...
              IG_A(4,1)+IG_A(2,3)])/2;
obs.Sm_B = ([IG_B(2,1)-IG_B(3,4), ...
             IG_B(3,1)-IG_B(4,2), ...
             IG_B(4,1)-IG_B(2,3)])/2;
obs.Sq_B = -([IG_B(2,1)+IG_B(3,4), ...
              IG_B(3,1)+IG_B(4,2), ...
              IG_B(4,1)+IG_B(2,3)])/2;

obs.F0 = F0;
obs.E0 = E0;
obs.S0 = S0;

end

function [C,t0] = build_block(chi,JJ0)
    % C: Wick's contractions of connected channel in
    % 3x3 spin average <sigma_i^{alpha} sigma_j^{alpha'}>_{\chi}.

    % t0: Simple Majorana MF hopping matrix t_{ij}^{\mu\nu} 
    % ( without gJ_{ij} and \partial gJ_{ij} / \partial \chi_{ij}^{\mu\nu} )

    c0 = chi(1,1);
    c1 = chi(1,2:4);
    c2 = chi(2:4,1);
    c3 = chi(2:4,2:4);
    u2 = [2 3 1]; u3 = [3 1 2];

    C = (c2*c1-c0*c3+ ...
         c3(:,u2).*c1(u3)-c3(:,u3).*c1(u2)+ ...
         c3(u2,:).*c2(u3)-c2(u2).*c3(u3,:)+ ...
         c3(u2,u3).*c3(u3,u2)-c3(u2,u2).*c3(u3,u3))/4;

    C1 = JJ0*c2; C2 = c1*JJ0; C3 = JJ0*c3; C4 = c3*JJ0;
    hh1 = C1+[C3(2,3)-C3(3,2);C3(3,1)-C3(1,3);C3(1,2)-C3(2,1)];
    hh2 = C2-[C4(2,3)-C4(3,2),C4(3,1)-C4(1,3),C4(1,2)-C4(2,1)];
    hh3_1 = -JJ0*c0;
    hh3_2 = JJ0(:,u3).*c1(u2)-JJ0(:,u2).*c1(u3)+ ...
            JJ0(u3,:).*c2(u2)-JJ0(u2,:).*c2(u3);
    hh3_3 = JJ0(u3,u2).*c3(u2,u3)+JJ0(u2,u3).*c3(u3,u2)- ...
            JJ0(u2,u2).*c3(u3,u3)-JJ0(u3,u3).*c3(u2,u2);
    t0 = [-sum(JJ0.*c3,'all'),hh1.';hh2.',hh3_1+hh3_2+hh3_3]/4;

end

function C = CoMatrix4(A)
% 4x4 cofactor matrix: C(i,j)=(-1)^(i+j)*det(A without i,j).
  C = zeros(4);

    for ii = 1:4
       for jj = 1:4     

           Bs = (-1)^(ii+jj);      
           B = A; B(ii,:) = []; B(:,jj) = [];      
           C(ii,jj) = Bs*det(B);  

       end       
    end
end