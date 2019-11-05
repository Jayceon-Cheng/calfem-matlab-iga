clc

clear all

tic
%%% 层合壳元
addpath(genpath('C:\Users\Administrator\Desktop\calfem-matlab-iga-master - (2)'));
count=1;


%  pro5;
B=cell(3,3);
B{1,1}=[-0.5,-0.5,0,1];
B{1,2}=[0,-0.5,0,1];
B{1,3}=[0.5,-0.5,0,1];
B{2,1}=[-0.5,0,0,1];
B{2,2}=[0,0,0,1];
B{2,3}=[0.5,0,0,1];
B{3,1}=[-0.5,0.5,0,1];
B{3,2}=[0,0.5,0,1];
B{3,3}=[0.5,0.5,0,1];
Xi=[0,0,0,1,1,1];
Eta=[0,0,0,1,1,1];
Zeta=[0,0,1,1];
p=numel(Xi)-size(B,1)-1;
q=numel(Eta)-size(B,2)-1;

deg.p = p; deg.q=q;
% deg.r = r;
clear p q

%% h-refinement (and solve the transform matrix for mmultigrid)

[ X_ ] = getSubDivKVValues( Xi, 4);
[B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,4);
[B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
% [ Z_ ] =getSubDivKVValues( Zeta, 1);
% [B,Zeta,Rz]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
KV.Xi=Xi; KV.Eta=Eta;

n1=size(B,1);m1=size(B,2);l1=size(B,3);
[INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );
B1=B(:,:,1:size(B,3));  %(9*3*12)  %无重合点时

% Number of control points after refinement
n = size(B1,1);
m = size(B1,2);
l = size(B1,3);
%   B1=B(:,:,1:l);
%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta;


% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) ; % number of elements

nnp = n*m; % number of global basis functions

nen = (deg.p+1)*(deg.q+1); % number of local basis functions
sdof = nnp*6; % number of global degrees of freedom
ldof = nen*6; % number of local degrees of freedom
% Build connectivity arrays
% [INN,IEN,AA] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
%   [INN,IEN,INB] = BLDINCIEN( deg,n,m,l );
[INN,IEN,INB] =BldINCIEN( deg,n,m,l );

ID = reshape(1:max(max(IEN))*6,6,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(6*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),6*nen,1);
end
ndof=size(IEN,1);  %单元节点数
%% Material parameters:
% thickness=0.1;   %plate thickness
density=1000;
nc=4 ;     %层数number of layers
h=0.24;     %总厚度
z= 0:h/nc:h;
zm = h/2;
zi = z-zm;
Zt = zeros(nc,2);
Zt(:,1)=zi(1:nc);
Zt(:,2)=zi(2:nc+1);

thetai=zeros(1,nc);
% thetai=ones(1,nc)*pi/180*90;

El = 137.9e9;       % 增大
Et = 10.34e9;       % 增大
vlt = 0.29;         % 减小
vtl = Et*vlt/El;
% vtl = 0.14;
Glt = 6.89e9;       % 增大
Glw = 6.89e9;
Gtw = 3.9e9;        % 增大
E=[El Et] ;   %沿纤维方向弹模1和横向弹模2
v=[vlt vtl];   %主泊松比vlt，
Q = [El/(1-vlt*vtl) vtl*El/(1-vlt*vtl) 0;
    vtl*El/(1-vlt*vtl) Et/(1-vlt*vtl)  0;
    0     0         Glt];
% %% Material parameters:
% E_Y = 210e4;
% nu_P = 0.3;
% damp=0.1;
% %D=hooke_strain(5,E_Y,nu_P);

%% Gauss-Legendre quadrature points:
[ gp_x,w_x ] = getGP( deg.p );
[ gp_y,w_y ] = getGP( deg.q );
% [ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
% NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation
Bmatrix=cell(size(IEN));
% Element loop
K = zeros(sdof); % Needs to be changed to sparse for large problems!!
M = zeros(sdof);
F = zeros(sdof,1);


for e =1:nel
    je=1;
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    %     nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)))
        continue
    end
    
    Ke = zeros(nen*6);
    Fe = zeros(nen*6,1);
    Me = zeros(nen*6);
    %     Ce = zeros(nen*3);
    
    %      for k = 1 : NQUADz
    %      for i = 1 : NQUADx % Loop trough Gauss points
    for j = 1 : NQUADy
        for i = 1 : NQUADx
            %               for k = 1 : NQUADz
            % Gauss point
            GP.xi_tilde = gp_x(i);
            GP.eta_tilde = gp_y(j);
            %GP.zeta_tilde = gp_z(k);
            
            % Get Basis, derivatives, and det(J) for current gauss pt
            [ R,dR_dx,Jdet ] =D2_Shape_function( GP,e,deg,B,KV,INN,IEN);
            % Combine quadrature weights with det(J)
            Jmod = abs(Jdet)*w_x(i)*w_y(j);
            
            % Build  Ke from membrane , dending ,shear components
            [ Ke_ ,B_m,dme] = shell_Build_K_Local( dR_dx,Jmod,R ,ndof,Zt,thetai,Q,nc,Gtw,Glw,h,density);
            %               [ Kg_ ,kinmtx ] = shell_Build_Kg_Local( dR_dx,Jmod,nen,thetai,R ,ndof);
            
            ns=[] ;
            for ii=1:ndof         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
                ns=[ns, 1+6*(ii-1):5+6*(ii-1)];
                Ke(6*ii,6*ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
            end
            Ke(ns,ns)= Ke(ns,ns)+Ke_;
            Me(ns,ns)= Me(ns,ns)+dme;
            Bm{je,e}=B_m; %membrane B matrix
            DR_dx{je,e}=dR_dx;
            je=je+1;
        end
    end
    %     end
    
    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
    %     KG(idx,idx) = KG(idx,idx) + Kg;
    M(idx,idx) = M(idx,idx) + Me;

end


%% Apply load in vertical direction on identified nodes that are listed in
% constNod:
constNod=[];
constNod=[constNod reshape(INB(:,1,:),1,numel(INB(:,1,:)))];
F(ID(1,constNod)) = 5e4;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
constNod = [];
constNod=[constNod reshape(INB(:,size(INB,2),:),1,numel(INB(:,size(INB,2),:)))];
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];

%% Solve system static
[a,K,M]=solveq(K,F,bc1,M);
% ex4(Xi,Eta,B,deg.p,deg.q )

%% calculate the contribuation of in-plane stress and construct the geometric matrix KG

ai=zeros(5*nen,nel);num_node=zeros(nnp,1);
stress_n=zeros(nnp,3);
Gb=zeros(2,5*ndof);
Gs1=zeros(2,5*ndof);
Gs2=zeros(2,5*ndof);
KGb=zeros(5*ndof);
KGs=zeros(5*ndof);
KG=zeros(sdof);
for sq = 1:nel
    ai(:,sq)=a(LM(ns,sq)) ;
    
    ni = INN(IEN(1,sq),1); nj = INN(IEN(1,sq),2);
    %       nk = INN(IEN(1,sq),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj))
        continue
    end
    gp=1;np=1;
    Kg=zeros(6*ndof);
    %     for k = 1 : NQUADz
    %         for i = 1 : NQUADx % Loop trough Gauss points
    for j = 1 : NQUADy
        for i = 1 : NQUADx
            
            kinmtx=Bm{gp,sq};
            strain(:,gp,sq)=kinmtx*ai(:,sq);
            siga=zeros(3,1);
            dN=DR_dx{gp,sq};
            Gb(1,1:5:5*ndof)=dN(:,1)';
            Gb(2,1:5:5*ndof)=dN(:,2)';
            Gs1(1,2:5:5*ndof)=dN(:,1)'; Gs1(2,2:5:5*ndof)=dN(:,2)';
            Gs2(1,3:5:5*ndof)=dN(:,1)'; Gs2(2,3:5:5*ndof)=dN(:,2)'; 
            for il = 1:nc    %layers
                thi = thetai(il);   %angles
                Ti=[cos(thi)^2 sin(thi)^2 -2*cos(thi)*sin(thi);
                    sin(thi)^2 cos(thi)^2 2*cos(thi)*sin(thi);
                    cos(thi)*sin(thi) -cos(thi)*sin(thi) cos(thi)^2-sin(thi)^2];
                z1 = Zt(il,1);   %Zt
                z2 = Zt(il,2);
                Qb=Ti*Q*Ti.';
                siga=siga+Qb*strain(:,gp,sq);
                Sigma=[siga(1),siga(3); siga(3),siga(2)];
                KGb=KGb+Gb'*Sigma*Gb*det(Jmod)*(z2-z1);
                KGs=KGs+(Gs1'*Sigma*Gs1+Gs2'*Sigma*Gs2)*det(Jmod)*(z2^3-z1^3)/3;
                
            end
            ns=[] ;
            for ii=1:ndof         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
                ns=[ns, 1+6*(ii-1):5+6*(ii-1)];
                %Kg(6*ii,6*ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
            end
            Kg(ns,ns)= Kg(ns,ns)+KGb+KGs;

            stress_n(IEN(gp,sq),:)=stress_n(IEN(gp,sq),:)+siga';
            %                  num_node(IEN(gp,sq))= num_node(IEN(gp,sq))+1;
            gp=gp+1;
        end
    end
    %     end
    idx = LM(:,sq)';
    KG(idx,idx) = KG(idx,idx) + Kg;
end

for ii=1:nel
    for jj=1:size(IEN,1)
        num_node(IEN(jj,ii))= num_node(IEN(jj,ii))+1;
    end
end
for i=1:n*m
    hh=INN(i,:);
    gcoord(i,:)=B1{hh(1),hh(2)}(1:3);
end
for i=1:nnp
    num=num_node(i);
    stress_n(i,:)=stress_n(i,:)./num;
end


[V,D]=eigs(K,KG,20,'sm');
[lamda,I]=sort(diag(D));
V=V(:,I);

toc

