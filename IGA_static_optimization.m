clc
clear all
tic

% li_IGA;
addpath(genpath('C:\Users\Administrator\Desktop\calfem-matlab-iga-master - (2)'));
%% 初始模型(直接)
main_IGA;  

%% 修改模型
ANSWER1 = questdlg('Select the example to run.', ...
                  'chose the type of calculation', ...
                  'time_analysis', 'static_reanalysis_local ', 'static_reanalysis_global', 'defualt');
%% 基于时态分析            
if strcmpi(ANSWER1,  'time_analysis')    
    %%----基于时态的程序
    return;
end
%% 基于重分析的静态结构优化
if strcmpi(ANSWER1,   'static_reanalysis_local ')

%% 基于重分析的优化 (local)
Choices = {'Algo1 (DE)', 'Algo2 (PSO)', 'Algo3 (SA)'};

ANSWER2 = questdlg('Select the example to run.', ...
                  'Growing Neural Gass', ...
                  Choices{1}, Choices{2}, Choices{3}, Choices{1});

if strcmpi(ANSWER2, Choices{1})
    de;
    return;
end

if strcmpi(ANSWER2, Choices{2})
    pso;
    return;
end

if strcmpi(ANSWER2, Choices{3})
    sa;
    return;
end
 return
end

if strcmpi(ANSWER1, 'static_reanalysis_global')
    my_HGWO_SVR;
    return;
end
%% 基于重分析的优化 (global)




