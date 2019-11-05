%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA102
% Project Title: Implementation of Particle Swarm Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% clc;
% clear;
% close all;
tic
load initial 
B1 = B0;
bB=cell2mat(B1);
B01=B1;

 nVar=(size(B1,2)-2)*size(B1,3);
lb=0.8*ones(1,nVar); % 参数取值下界
ub=1.2*ones(1,nVar);% 参数取值上界  

[INBz1,INNz1,Br1,Bz11,Bz] =Midgridx(B0,lb,ub,Max,Min,KV,deg) ;
[ Pp1,nsub,ssub,idsub]=suoyin_1(INBz1,INNz1,Br1,Bz11,B01,Bz); %1初始,2参数化修改
VarMax=ub;
VarMin=lb;
% %% Problem Definition
% 
% CostFunction=@(x) REiga_pro3(x);        % Cost Function
% 
% % nVar=1;            % Number of Decision Variables
% % 
% % VarSize=[1 nVar];   % Size of Decision Variables Matrix
% 
% nVar1 = 1;    % Number of Decision Variables
% 
% % nVar2 = 1;
%  VarSize2 = [1 nVar1];
% VarSize1 = [1 nVar1];
% 
% % VarMin=1.15;        % Lower Bound of Variables
% % VarMax=1.4;         % Upper Bound of Variables
% 
% V1Min = 1.18;   % Lower Bound of Variables
% V1Max = 1.26;    % Upper Bound of Variables
% V2Min = 1.9;   % Lower Bound of Variables
% V2Max = 2.3;

%% PSO Parameters

MaxIt=4;      % Maximum Number of Iterations

nPop=4;        % Population Size (Swarm Size)

% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1.*(ub-lb);
VelMin=-VelMax;
% VelMax1=0.1*(VarMax1-VarMin1);
% VelMin1=-VelMax1;
% VelMax2=0.1*(VarMax2-VarMin2);
% VelMin2=-VelMax2;
% VelMin=[VelMin1,VelMin2];
% VelMax=[VelMax1,VelMax2];
%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
%     x=init_individual(lb,ub,nVar,nPop); % 随机初始化位
%     x1=unifrnd(VarMin1,VarMax1,VarSize1);
%     x2=unifrnd(VarMin2,VarMax2,VarSize2);\
xRange=ub-lb;
xLower=lb;
x=rand(1,nVar).*xRange+xLower;
%     
%     x=[x1,x2];
    particle(i).Position=x;
    % Initialize Velocity
%     v1=zeros(VarSize1);v2=zeros(VarSize2);v=[v1,v2];
    v=zeros(1,nVar);
    particle(i).Velocity=v;
    % Evaluation
%     particle(i).Cost=CostFunction(particle(i).Position);
   particle(i).Cost=objecfun(particle(i).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(1,nVar).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(1,nVar).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
%         particle(i).Cost = CostFunction(particle(i).Position);
        particle(i).Cost =objecfun(particle(i).Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
%         % Plot F1 Costs
%     figure(1);
%     PlotCosts(GlobalBest,it);
%     pause(0.01);
%     hold on
    w=w*wdamp;
    
end
toc
BestSol = GlobalBest;

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;

