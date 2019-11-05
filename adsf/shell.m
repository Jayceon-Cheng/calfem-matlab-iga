function [Ke,Me] = shell(i)
global thetai nodes gcoord coord ang nl
p = gcoord(nodes(i,:),:);
dp(1,:) = p(2,:)-p(1,:);
dp(2,:) = p(3,:)-p(1,:);
dp(3,:) = p(4,:)-p(1,:);
x = dp(1,:)/sqrt(dp(1,:)*(dp(1,:))');
z1 = cross(dp(1,:),dp(2,:));
z = z1/sqrt(z1*z1');
y = cross(z,x);
coord = zeros(4,2);
local = [x',y'];
coord(2:4,:) = dp*local;
% angi = ang(i,:);
% for i = 1:nl
%     t = angi(i);
%     xf = [cos(t) sin(t) 0];
%     zf = cross(x,xf);
%     thetai(i) = sign(zf*z')*acos(xf*x');
% end
thetai = ang(i,1:nl);
[ke1,me1] = elematoflaminate;
ke = zeros(24,24);
me = zeros(24,24);
a = 1:6:19;
aax = bsxfun(@plus,a,(0:4)')';
bbx = bsxfun(@plus,a,(0:4)');
aax = aax(:)';
bbx = bbx(:)';
ke(aax,aax) = ke(aax,aax)+ke1;
me(bbx,bbx) = me(bbx,bbx)+me1;

for ii=6:6:24         % 一个节点所在多个单元处于同一平面内时，刚度矩阵会奇异，这里给第六自由度(对角)赋值消除奇异
  ke(ii,ii) = 0.05;  % 圆柱壳在0.05附近,越小频率越小,小于0.01不稳定；板对大小不敏感
end

% te=zeros(6);
% te(1:3,4:6)=eye(3);
% te(4:6,1:3)=eye(3);
% % te(6,6)=1;
% Te=blkdiag(te,te,te,te);
% ke=Te*ke*Te';
% me=Te*me*Te';
te=zeros(6);
te(1,4)=1;
te(2,5)=1;
te(3,1)=1;
te(4,3)=1;
te(5,2)=1;
te(6,6)=1;
Te=blkdiag(te,te,te,te);
ke=Te*ke*Te';
me=Te*me*Te';

% X = [1 0 0];
% Y = [0 1 0];
% Z = [0 0 1];
% lx = x*X';
% mx = x*Y';
% nx = x*Z';
% ly = y*X';
% my = y*Y';
% ny = y*Z';
% lz = z*X';
% mz = z*Y';
% nz = z*Z';
% T3 = [lx mx nx;
%       ly my ny;
%       lz mz nz];
T3 = [x;y;z];
T = blkdiag(T3,T3,T3,T3,T3,T3,T3,T3);
Ke = T'*ke*T;
Me = T'*me*T;

