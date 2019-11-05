% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

clc;
clear;
close all;

load initial 

bB=cell2mat(B1);
B01=B1;

 nVar=size(B1,2)*size(B1,3);
lb=0.8*ones(1,nVar); % 参数取值下界
ub=1.2*ones(1,nVar);% 参数取值上界  
[INBz1,INNz1,Br1,Bz11,Bz] =Midgridx(B0,lb,ub,Max,Min,KV,deg) ;
[ Pp1,nsub,ssub,idsub]=suoyin_1(INBz1,INNz1,Br1,Bz11,B01,Bz); %1初始,2参数化修改
VarMax=ub;
VarMin=lb;

%% Problem Definition
tic
% CostFunction = @(x) REiga_pro3(x);  % Cost Function
CostFunction = @(x)objecfun(x,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
% nVar = 1;    % Number of Decision Variables
%  nVar2 = 1;
%  nVar3 = 2;

VarSize1 = [1 nVar]; % Size of Decision Variables Matrix
 VarSize2 = [1 nVar];
% % VarSize3 = [1 nVar3];

% V1Min = 1.15;   % Lower Bound of Variables    %benchmark
% V1Max = 1.4;    % Upper Bound of Variables
% V1Min = 1.18;   % Lower Bound of Variables    %screwdriver
% V1Max = 1.26;    % Upper Bound of Variables
% V2Min = 1.9;   % Lower Bound of Variables    %screwdriver
% V2Max = 2.3;    % Upper Bound of Variables

% x1=unifrnd(V1Min, V1Max, VarSize1);
% x2=unifrnd(V2Min, V2Max, VarSize2);
% x=[x1,x2];
xRange=ub-lb;
xLower=lb;
x=rand(1,nVar).*xRange+xLower;
% Number of Objective Functions
%nObj = numel(CostFunction(unifrnd(V1Min, V1Max, VarSize1)));  %ladingcouyang
nObj = numel(CostFunction(x));

%% NSGA-II Parameters

% Generating Reference Points
nDivision = 10;
Zr = GenerateReferencePoints(nObj, nDivision);

MaxIt = 1;  % Maximum Number of Iterations

nPop = 2;  % Population Size

pCrossover = 0.8;       % Crossover Percentage
nCrossover = 2*round(pCrossover*nPop/2); % Number of Parnets (Offsprings)

pMutation = 0.3;       % Mutation Percentage
nMutation = round(pMutation*nPop);  % Number of Mutants

mu = 0.02;     % Mutation Rate
VarMax=ub;
VarMin=lb;
sigma = 0.1*(VarMax-VarMin); % Mutation Step Size
% sigma = 0.1*(V1Max-V1Min);

%% select Parameters

params.nPop = nPop;
params.Zr = Zr;
params.nZr = size(Zr,2);
params.zmin = [];
params.zmax = [];
params.smin = [];

%% Initialization

disp('Staring NSGA-III ...');

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.NormalizedCost = [];
empty_individual.AssociatedRef = [];
empty_individual.DistanceToAssociatedRef = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    %pop(i).Position = unifrnd(VarMin, VarMax, VarSize); % varibles
%      x1=unifrnd(V1Min, V1Max, VarSize1);
%      x2=unifrnd(V2Min, V2Max, VarSize2);
% % %      x3=unifrnd(V3Min, V3Max, VarSize3);
%       x=[x1 x2];
      xRange=ub-lb;
xLower=lb;
x=rand(1,nVar).*xRange+xLower;
%  x=x1;
    pop(i).Position = x;
    pop(i).Cost = CostFunction(pop(i).Position);
end

% Sort Population and Perform Selection
[pop, F, params] = SortAndSelectPopulation(pop, params);

disp('Initialization.');
%% NSGA-II Main Loop

for it = 1:MaxIt
 
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2

        i1 = randi([1 nPop]);
        p1 = pop(i1);

        i2 = randi([1 nPop]);
        p2 = pop(i2);

        [popc(k, 1).Position, popc(k, 2).Position] = Crossover(p1.Position, p2.Position,VarMax,VarMin);

        popc(k, 1).Cost = CostFunction(popc(k, 1).Position);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position);

    end
    popc = popc(:);

    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation

        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, mu, sigma,VarMax,VarMin);

        popm(k).Cost = CostFunction(popm(k).Position);

    end
 
    % Merge
    pop = [pop
           popc
           popm]; %#ok
    
    % Sort Population and Perform Selection
    [pop, F, params] = SortAndSelectPopulation(pop, params);
    
    % Store F1
    if numel(F{1})==1
    F1 = pop(F{1});
    else
     F1 = pop(F{1}(1)); 
    end
   
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Cost = ' num2str(F1.Cost)]);

    % Plot F1 Costs
    figure(1);
    PlotCosts(F1,it);
    pause(0.01);
    hold on
    GA(it,1:nVar)=F1.Position;
    GA(it,nVar+1)=F1.Cost;
    Error(it)=error_objecfun(F1.Position,nVar,INBz1,INNz1,Br1,Bz11,Bz,nsub,ssub,idsub,Pp1);
   
end
% reana(F1);

%% Results

disp(['Final Iteration: Number of F1 Members = ' num2str(numel(F1))]);
disp('Optimization Terminated.');
toc

