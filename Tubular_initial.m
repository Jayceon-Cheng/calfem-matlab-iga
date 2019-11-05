
clear all
clc
clf
tic
%% Load data
addpath(genpath('..\calfem-matlab-iga-master-screwdriver'));
count=1;

%  data_l;
%      run ./../NURBS/Solid/data_solid_chilun
%           run ./../NURBS/solid/data_solid9
%           run ./../NURBS/Solid/data_part2
%             run ./../NURBS/Solid/data_multi_part
      
%  for i=1:size(B,1)
%   for j=1:size(B,2)
%         B{i,j,3}(2)=B{i,j,3}(2)+0.5;
%   end
% end   
%}

 pro5;
 %轴向坐标不变，将内部圆柱面变为均匀一致性；
 for i=1:size(B,2)
     
      B{2,i,1}(2:3)=[1;0];
      B{2,i,2}(2:3)=[1;1];
      B{2,i,3}(2:3)=[0;1];
      B{2,i,4}(2:3)=[-1;1];
      B{2,i,5}(2:3)=[-1;0];
      B{2,i,6}(2:3)=[-1;-1];
      B{2,i,7}(2:3)=[0;-1];
      B{2,i,8}(2:3)=[1;-1];
      B{2,i,9}(2:3)=[1;0];
 end

% load Pso BestSol
% for i=1:numel(B)
%     B{i}(1:2:3)=B{i}(1:2:3)./10;
% end
% B=permute(B,[3 2 1]);
deg.p = p; deg.q=q; deg.r = r; clear p q r

%% h-refinement (and solve the transform matrix for mmultigrid) 

[ X_ ] = getSubDivKVValues( Xi, 1); 
[B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,1);
[B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta,Rz]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;

i=1:size(Rx,1);j=1:size(Ry,1); k=1:size(Rz,1);
s=1:size(Rx,2);t=1:size(Ry,2);l=1:size(Rz,2);
[ii,jj,kk]=meshgrid(i,j,k); [ss,tt,ll]=meshgrid(s,t,l);

cn=ii+(jj-ones(size(jj)))*size(Rx,1)+(kk-ones(size(jj)))*size(Ry,1)*size(Rx,1); 
xn=ss+(tt-ones(size(tt)))*size(Rx,2)+(ll-ones(size(ll)))*size(Ry,2)*size(Rx,2);
zz=kron(Rx(i,s,size(Rx,3)),Ry(j,t,size(Ry,3)));
RR=kron(zz,Rz(:,:,size(Rz,3)));
 RR=kron(RR,eye(3));  %疏密网格间的限制矩阵
 Rs=sparse(RR);
 P1=RR';          %疏密网格间的插值矩阵
clear RR cn xn zz ii jj kk ss tt ll
numel(Rs(Rs~=0));

%  [ X_ ] = getSubDivKVValues( Xi, 1); 
%  [B,Xi,Rx]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
%   [ E_ ] = getSubDivKVValues( Eta,2);
%   [B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta,Rz]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
  Rx=eye(size(Rx,2));
   Ry=eye(size(Ry,2));
KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;
%对疏密网格自由度进行全局编号《得到总体自由度的转换矩阵的元素为
i=1:size(Rx,1);j=1:size(Ry,1);k=1:size(Rz,1);
s=1:size(Rx,2);t=1:size(Ry,2);l=1:size(Rz,2);
[ii,jj,kk]=meshgrid(i,j,k);
[ss,tt,ll]=meshgrid(s,t,l);
% [rx,ry,rz]=meshgrid(Rx(:,:,size(Rx,3)),Ry(:,:,size(Ry,3)),Rz(:,:,size(Rz,3)));

cn=ii+(jj-ones(size(jj)))*size(Rx,1)+(kk-ones(size(jj)))*size(Ry,1)*size(Rx,1); 
xn=ss+(tt-ones(size(tt)))*size(Rx,2)+(ll-ones(size(ll)))*size(Ry,2)*size(Rx,2);
zz=kron(Rx(i,s,size(Rx,3)),Ry(j,t,size(Ry,3)));
RR=kron(zz,Rz(:,2:size(Rz,2),size(Rz,3)));
 RR=kron(RR,eye(3));  %疏密网格间的限制矩阵
 Rs2=sparse(RR);
P2=RR';          %疏密网格间的插值矩阵
clear RR cn xn zz ii jj kk ss tt ll
numel(Rs(Rs~=0));



%{
% Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
[ X_ ] = getSubDivKVValues( Xi, 1); 
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta,2);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 2);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );
%}

KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;
BJ=cell2mat(B);
for i=1:size(B,1)
    for j=1:size(B,2)
     Norm(i,j)=   norm(B{i,j,1}(2:3));
    end
end
Max=round(max(Norm(1,:)));
if Max<max(Norm(1,:))
 Max=Max+1;
end
Min=round(min(Norm(3,:)));
if Min>min(Norm(3,:))
    Min=Min-1;
end
%  Xmax=max(max(max(BJ(1:4:end,:,:))))+1;
%  Xmin=min(min(min(BJ(1:4:end,:,:))))-1;
%  Ymax=max(max(max(BJ(2:4:end,:,:))))+1;
%  Ymin=min(min(min(BJ(2:4:end,:,:))));
%  Zmax=max(max(max(BJ(3:4:end,:,:))))+1;
%  Zmin=min(min(min(BJ(3:4:end,:,:))));
   n1=size(B,1);m1=size(B,2);l1=size(B,3);
  [INN1,IEN1] = BldINCIEN( deg,n1,m1,l1 );
B1=B(:,:,2:size(B,3));  %(9*3*12)
%   nVar=180;
%  parent_Position1=1.2*ones(1,180);
% % 
% % parent_Position1=PA;
% parent_Position1=1.2*ones(1,180);
% parent_Position1(180-11)=0.8;
% for i=1:180
%     B1{1+3*(i-1)}(2:3)=parent_Position1(i).*B1{1+3*(i-1)}(2:3);
%     B1{2+3*(i-1)}(2:3)=parent_Position1(i).*B1{2+3*(i-1)}(2:3);
% %     B1{3*(i)}(2:3)=parent_Position1(i).*B1{3*(i)}(2:3);
% 
% % end
% load posT
%     parent_Position1=post(495,:);
% 
% for i=1:180
%     B1{1+3*(i-1)}(2:3)=parent_Position1(i).*B1{1+3*(i-1)}(2:3);
%     B1{2+3*(i-1)}(2:3)=parent_Position1(i).*B1{2+3*(i-1)}(2:3);
% %     B1{3*(i)}(2:3)=parent_Position1(i).*B1{3*(i)}(2:3);
% 
% end
% B(:,:,2:size(B,3))=B1;
% B(:,:,1)=B(:,:,size(B,3));
% [ E_ ] = getSubDivKVValues( Eta,2);
% [B,Eta,Ry]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
% KV.Xi=Xi; KV.Eta=Eta; KV.Zeta=Zeta;
% B1=B(:,:,2:size(B,3));  %(9*3*12)

% Number of control points after refinement
n = size(B1,1);
m = size(B1,2);
l = size(B1,3);
%   B1=B(:,:,1:l);
%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) * (l1-deg.r); % number of elements

nnp = n*m*l; % number of global basis functions

nen = (deg.p+1)*(deg.q+1)*(deg.r+1); % number of local basis functions
ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
% [INN,IEN,AA] = BldINCIEN( deg,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
 [INN,IEN,INB] = BLDINCIEN11( deg,n,m,l );
% for i=1:numel(IEN)
%     ie=IEN(i);
%    if (ie>=(n*m1+1)&&ie<=n*m)||(ie>=(n*(2*m-1)+1)&&ie<=n*m*2)||(ie>=(n*(l*m-1)+1)&&ie<=n*m*l)
%           IEN(i)=IEN(i-n*m1);
%    end
%   
% end
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed) 在系统中单元节点或者说是自由度编号。
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end

%% Material parameters:
E_Y = 210e4;
nu_P = 0.3;
denisty=785;
damp=0.1;
D=hooke_strain(5,E_Y,nu_P);

%% Gauss-Legendre quadrature points: 
[ gp_x,w_x ] = getGP( deg.p );     
[ gp_y,w_y ] = getGP( deg.q );
[ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation
Bmatrix=cell(size(IEN));
% Element loop
K = zeros(ndof); % Needs to be changed to sparse for large problems!!
% K = zeros(ndof);
F = zeros(ndof,1);
for e =1:nel
    je=1;
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk)))&&e<=(n-deg.p) * (m-deg.q) * (l-deg.r)
        continue
    end
   
    
    Ke = zeros(nen*3);
    Fe = zeros(nen*3,1);
%     Me = zeros(nen*3);
%     Ce = zeros(nen*3);
%     
%     for i = 1 : NQUADx % Loop trough Gauss points
%         for j = 1 : NQUADy
%             for k = 1 : NQUADz
     for k = 1 : NQUADz
%      for i = 1 : NQUADx % Loop trough Gauss points
         for j = 1 : NQUADy
             for i = 1 : NQUADx 
%               for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,e,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_,kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                Ke = Ke + Ke_;
                Bmatrix{je,e}=kinmtx; 
                je=je+1;
%                 % Build Me
%                  Ni=eye(3);
%                  for s=1:length(R)
%                      shamtx(:,(s-1)*3+1:s*3)=R(s)*Ni;
%                  end
%                  Me_=denisty*shamtx'*shamtx*Jmod;
%                  Ce_=damp*shamtx'*shamtx*Jmod;
%                  Me = Me + Me_;
%                  Ce = Me + Me_;
            end
        end
    end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
%     M(idx,idx) = M(idx,idx) + Me;
%     C(idx,idx) = C(idx,idx) + Ce;
%     F(idx) = F(idx)+Fe;
end
% Apply load at two corner nodes
%     loc_z=0;loc_x=0;   %data_10的
%     loc_y=0;loc_x=0;
%   loc_y=8;loc_x=3;  %chilun
% loc_x=0; loc_y=0.0064; %guan
loc_x=0; loc_y=0; %guan
constNod1 = [];
% for i = 1 : numel(B1)
%          if B1{i}(2) == loc_y||abs(B1{i}(2)-0.0064)<=1e-4
% %            if (B{i}(2) == 0.495||B{i}(2) == 0.505)
%              if B1{i}(1) == loc_x||abs(B1{i}(1)-0)<=1e-4
% %              if (B1{i}(1) == loc_x||abs(B1{i}(1))<=1e-10)
%                 constNod=[constNod i];
%              end
%           end
%         
% %         end
% end

% Apply load in vertical direction on identified nodes that are listed in
% constNod:
constNod1=[];
%  constNod=[constNod reshape(INB(:,1,l/4+1),1,numel(INB(:,1,l/4+1)))];
%  constNod=[constNod reshape(INB(:,1,3*l/4+1),1,numel(INB(:,1,3*l/4+1)))];

 constNod1=[constNod1 reshape(INB(1,:,size(INB,3)),1,numel(INB(1,:,size(INB,3))))];
%  constNod1=[constNod1 reshape(INB(:,1,1:19),1,numel(INB(:,1,1:19)))];
F(ID(2,constNod1)) = -2e5;
% constNod2=[];
% constNod2=[constNod2 reshape(INB(:,1,20),1,numel(INB(:,1,20)))];
% F(ID(1,constNod2)) =1.5e4;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
%   loc_z=-1.2;
%   loc_y=-1.3;
   loc_y=20;
constNod = [];
% for i = 1 : numel(B1)
%     if B1{i}(2) == loc_y
%         constNod=[constNod i];
%     end
% end
constNod=[constNod reshape(INB(:,size(INB,2),:),1,numel(INB(:,size(INB,2),:)))];
constNod=[constNod reshape(INB(:,1,:),1,numel(INB(:,1,:)))];
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];


%% Solve system
[a,K,F]=solveq(K,F,bc);

%%
ai=zeros(3*nen,nel);num_node=zeros(nnp,1);
stress_n=zeros(nnp,4);
for sq = 1:nel
     ai(:,sq)=a(LM(:,sq)) ;
    
      ni = INN(IEN(1,sq),1); nj = INN(IEN(1,sq),2);nk = INN(IEN(1,sq),3);
    
    % Check if element has zero measure
    if ((KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk)))&&sq<=(n-deg.p) * (m-deg.q) * (l-deg.r)
        continue
    end
    gp=1;np=1;
%     for i = 1 : NQUADx % Loop trough Gauss points
%         for j = 1 : NQUADy
%             for k = 1 : NQUADz
    for k = 1 : NQUADz
%         for i = 1 : NQUADx % Loop trough Gauss points
         for j = 1 : NQUADy
             for i = 1 : NQUADx 
%                 % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ R,dR_dx,Jdet ] = Shape_function( GP,sq,deg,B1,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                [ Ke_, kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen );
                
%                 kinmtx=Bmatrix{gp,sq};
                strain(:,gp,sq)=kinmtx*ai(:,sq);
                siga=D*strain(:,gp,sq);
                
                stress(gp,sq,1:3)=siga(1:3);
          % von Mises stress
                stress(gp,sq,4)  = sqrt(((siga(1)-siga(2))^2+(siga(3)-siga(2))^2+...
                                   (siga(1)-siga(3))^2)/2);
                               aaa(1,:)=stress(gp,sq,:);  
                stress_n(IEN(gp,sq),:)=stress_n(IEN(gp,sq),:)+aaa(1,:);
%                 num_node(IEN(gp,sq))= num_node(IEN(gp,sq))+1;
                gp=gp+1;
            end
        end
    end   
end

for ii=1:nel
    for jj=1:size(IEN,1)
        num_node(IEN(jj,ii))= num_node(IEN(jj,ii))+1;
    end
 end
for i=1:n*m*l
  h=INN(i,:)  ;
gcoord(i,:)=B1{h(1),h(2),h(3)}(1:3);
end
for i=1:nnp
    num=num_node(i);
    stress_n(i,:)=stress_n(i,:)./num;
end

%%------------环形模型需增加的的应力操作---------
%{
 mise=zeros(n1*m1*l1,4);ip=size(B1,2)*size(B1,1);
 for i=1:ip
     mise((l1*(i-1)+1):(l1*i-1),:)=stress_n(((l1-1)*(i-1)+1):(l1-1)*i,:);
     mise(l1*i,:)=stress_n((l1-1)*(i-1)+1,:);
 end
%}
%  mise=zeros(n1*m1*l1,4);ip=size(B1,2)*size(B1,1);
%  
%  mise(1:n1*m1*(l1-1),:)=stress_n;
%  mise((n1*m1*(l1-1)+1):end,:)=stress_n(1:ip,:);
 

%% TRANSFORM the matrix of element(p+1)*(q+1)(r+1) into 2*2*2
 e=0;bi=0;
    p= deg.p ; q=deg.q; r=deg.r ;
for k = 1 : (l)
        for j = 1 : (m-1)
            for i = 1 : (n-1)
                e=e+1;
               for kloc = 0 : 1   %1 ,  3
                   for jloc = 0 : 1   %2 ,  12
                           for iloc = 0 : 1    %3 ,  55                     
                              if (kloc+k)==(l+1)
                               
                                B_j=INB((iloc+i), jloc+j,1); 
                               else
                                B_j=INB((iloc+i), jloc+j, kloc+k);   % global function number
                              end
                              
                              b = (1+1)*(1+1)*(1+1)-kloc*(1+1)*(1+1) - jloc*(1+1)- iloc ; % local function number
                              BR(b,e) = B_j; % assign connectivity
                           end
                   end
               end
             
            end
        end
end
    BR1=BR;
    BR1(1,:)=BR(2,:);
    BR1(2,:)=BR(1,:);
    BR1(5,:)=BR(6,:);
    BR1(6,:)=BR(5,:);
    
%   tecplot(BR1,gcoord,mise);
     tecplot(BR1,gcoord,stress_n,a);%     tecplot(BR1,gcoord,stress_n);



%%
u = cell(size(B1)); 
comb = cell(size(B1));
for i = 1 : size(ID,2) %
u{i} = [a(ID(:,i)); 0]; 
comb{i} = B1{i} + u{i}; 
end
       comb(:,:,l1)=comb(:,:,1);
       comb1=cell2mat(comb);
%      tecplot(BR1,gcoord,stress_n,a);%     tecplot(BR1,gcoord,stress_n);
%% plot

plotNurbsSolidElementSimple( KV, B )
title('Initial geometry_the circular plate coupled to a doubly curved barrel shape')   %data_solid10

 plotNurbsSolidElementSimple( KV, comb )
 title('Displacements')
% 
 %}
 disp0=a;K0=K; ndof0=ndof;B0=B;

  save initial  disp0 ndof0  K0  B0 deg KV D B1 P1 P2 Max Min D LM  ID  INN IEN INB F stress_n
toc

