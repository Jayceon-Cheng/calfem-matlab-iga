function [  ] = plotNurbsCurve2DSimple( KV, B ) 
%[  ] = plotNurbsCurve2DSimple( curve )
%-------------------------------------------------------------
% PURPOSE:
%   Plot NURBS curve geometry 2d, and draw knots (elements).
%
%
% INPUT: Curve = Curve with:
%        Curve.KV (knot vector)
%        Curve.w (column vector storing weights)
%        Curve.CP (Matrix where each row contains a point)
%
% OUTPUT: none
%-------------------------------------------------------------



%% Generate basis

res = 100; % Compute 40 samples, uniformly distributed (and knot points are auto-included)
% [ R , U] = nrbasis_num( KV.res = 100; % Compute 40 samples, uniformly distributed (and knot points are auto-included)
[ R , U] = nrbasis_num( KV.Xi, B, res );

%% Calculate NURBS-spline curve
% C = R' * Curve.CP;
for i = 1:size(R,1)
    for j = 1 : size(R,2)
        
            Cx = Cx + R{i,j,k}*B{i,j,k}(1);
            Cy = Cy + R{i,j,k}*B{i,j,k}(2);

       
    end
end

%% Plot basis
%figure(1)
%plotNurbsBasis(R,U);

%% PLot curve
%figure(2)
plotNurbsCurve2D(C(:,1),C(:,2),Curve.KV,U,Curve.CP(:,1),Curve.CP(:,2));

end

