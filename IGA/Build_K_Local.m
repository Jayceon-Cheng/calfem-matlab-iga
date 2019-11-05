function [Ke, kinmtx ] = Build_K_Local( dR_dx,Jmod,D,nen )
% [ Ke ] = Build_K_Local( dR_dx,Jmod,D,nen )
%-------------------------------------------------------------
% PURPOSE:
% Calculate the contribution from the current integration point 
% to the linear elastic stiffness matrix for a solid 
% NURBS element.
%
% INPUT: dR_dx = vector of basis function derivatives (nen x 3)
%
%        Jmod  = Jacobian determinant
%
%        D     = constitutive matrix (6 x 6)
%
%        nen   = number of local basis functions
%
% OUTPUT: Ke = element stiffness matrix contribution (nen*3 x nen*3)
%-------------------------------------------------------------

% Generate B matrix
kinmtx=zeros(6,nen*3);
kinmtx(1,1:3:end) = dR_dx(:,1);
kinmtx(2,2:3:end) = dR_dx(:,2);
kinmtx(3,3:3:end) = dR_dx(:,3);
% kinmtx(4,2:3:end) = dR_dx(:,3);
% kinmtx(4,3:3:end) = dR_dx(:,2);
kinmtx(4,1:3:end) = dR_dx(:,2);
kinmtx(4,2:3:end) = dR_dx(:,1);
% kinmtx(5,1:3:end) = dR_dx(:,3);
% kinmtx(5,3:3:end) = dR_dx(:,1);
kinmtx(5,2:3:end) = dR_dx(:,3);
kinmtx(5,3:3:end) = dR_dx(:,2);
% kinmtx(6,1:3:end) = dR_dx(:,2);
% kinmtx(6,2:3:end) = dR_dx(:,1);
kinmtx(6,1:3:end) = dR_dx(:,3);
kinmtx(6,3:3:end) = dR_dx(:,1);
%stiffness matrix
Ke = kinmtx'*D*kinmtx*Jmod;
end

