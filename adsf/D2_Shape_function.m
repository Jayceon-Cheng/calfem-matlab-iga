function [ R,dR_dx,J]= D2_Shape_function( GP,e,deg,B,KV,INN,IEN)


p = deg.p ;   q = deg.q  ; 

% number of local basis functions:
nen = (p+1)*(q+1);
 %  nen = (p+1)*(q+1);                       
% NURBS coordinates; convention consistent with Algorithm 7
ni = INN(IEN(1,e),1);
nj = INN(IEN(1,e),2);
% if e>(size(IEN,2)-size(B,1)*size(B,2))
%  nk = size(B,3);   
% else
% nk = INC(IEN(1,e),3);
% end

% Calculate parametric coordinates from parent element coordinates
% Knot vectors KV.Xi, KV.Eta, and KV.Zeta and
% parent element coordinates xi_tilde, eta_tilde, zeta_tilde
% are given as input
  
xi = ((KV.Xi(ni+1)-KV.Xi(ni))*GP.xi_tilde ...
+ (KV.Xi(ni+1)+KV.Xi(ni))) / 2;
eta = ((KV.Eta(nj+1)-KV.Eta(nj))*GP.eta_tilde ...
+ (KV.Eta(nj+1)+KV.Eta(nj))) / 2;
% zeta = ((KV.Zeta(nk+1)-KV.Zeta(nk))*GP.zeta_tilde ...
% + (KV.Zeta(nk+1)+KV.Zeta(nk))) / 2;

% Calculate univariate B-spline functions using (2.1) and (2.2)
% and their derivatives using (2.12)
N1 = Der1BasisFun(ni-1,xi,p,KV.Xi)'; % xi-dir.
N2 = Der1BasisFun(nj-1,eta,q,KV.Eta)'; % eta-dir.
% N3 = Der1BasisFun(nk-1,zeta,r,KV.Zeta)'; % zeta-dir.
N = N1(:,1);     %三维的形函数N,M,L
dN_dxi = N1(:,2);
M = N2(:,1);
dM_deta = N2(:,2);
% L = N3(:,1);
% dL_dzeta = N3(:,2);
clear N1 N2 

% Build numerators and denominators (in local numbering)
x = zeros(1,nen); y = zeros(1,nen);
% z = zeros(1,nen);
R = zeros(nen,1);   % Array of trivariate NURBS basis functions
dR_dxi = zeros(nen,2); % Trivariate NURBS function derivatives
                       % w.r.t. parametric coordinates
loc_num = 0; % Local basis function counter
%  for k = 0 : r
    for j = 0 : q
        for i = 0 : p
            loc_num = loc_num + 1;
         
           R(loc_num) = N(p+1-i)*M(q+1-j) ...
                        * B{ni-i,nj-j}(4); % Function numerator (N*M*L*w)
            
            x(loc_num) = B{ni-i,nj-j}(1) ;
            y(loc_num) = B{ni-i,nj-j}(2);
            z(loc_num) = B{ni-i,nj-j}(3);
          
            %w(loc_num) = B{ni-i,nj-j,nk-k}(4);
            
           dR_dxi(loc_num,1) = dN_dxi(p+1-i)*M(q+1-j)...
           * B{ni-i,nj-j}(4); % Derivative numerator (dN*M*L*w)
%              dR_dxi(loc_num,1) = dN_dxi(p+1-i)*M(q+1-j) ...
%             * B{ni-i,nj-j}(4); % Derivative numerator (dN*M*L*w)
        
            dR_dxi(loc_num,2) = N(p+1-i)*dM_deta(q+1-j) ...
            * B{ni-i,nj-j}(4); % Derivative numerator (N*dM*L*w)
% 
%            dR_dxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dL_dzeta(r+1-k) ...
%            * B{ni-i,nj-j,nk-k}(4); % Derivative numerator (N*M*dL*w)
        end
        
     end
%  end
W = sum(R); % Function denominator (Sum(N*M*L*w))
dW_dxi = sum(dR_dxi(:,1)); % Derivative denominator (Sum(dN*M*L*w))
dW_deta = sum(dR_dxi(:,2)); % Derivative denominator (Sum(N*dM*L*w))
% dW_dzeta = sum(dR_dxi(:,3)); % Derivative denominator (Sum(N*M*dL*w))
            
% Divide by denominators to complete definitions of functions
% and derivatives w.r.t. parametric coordinates
dR_dxi(:,1) = (dR_dxi(:,1)*W - R*dW_dxi) / W^2;         
dR_dxi(:,2) = (dR_dxi(:,2)*W - R*dW_deta) / W^2;
% dR_dxi(:,3) = (dR_dxi(:,3)*W - R*dW_dzeta) / W^2;
R = R/W;

% Gradient of mapping from parameter space to physical space
% dx_dxi = [x;y;z] * dR_dxi;
dx_dxi = [x;y] * dR_dxi;

% Compute derivatives of basis functions
% with respect to physical coordinates
dR_dx = dR_dxi/dx_dxi; %dR_dxi * inv(dx_dxi)

% Gradient of mapping from parent element to parameter space
dxi_dtildexi=zeros(2); % Derivative of parametric coordinates
                       % w.r.t. parent element coordinates
dxi_dtildexi(1,1) = (KV.Xi(ni+1)-KV.Xi(ni))/2;
dxi_dtildexi(2,2) = (KV.Eta(nj+1)-KV.Eta(nj))/2;
% dxi_dtildexi(3,3) = (KV.Zeta(nk+1)-KV.Zeta(nk))/2;

% Compute the Jacobian
J_mat = dx_dxi*dxi_dtildexi;  %(3*3)  %取决于模型维度，size=（维度，维度）;

% Compute Jacobian determinant
J = det(J_mat);   

end

