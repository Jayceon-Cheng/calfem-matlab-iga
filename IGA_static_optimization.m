clc
clear all
tic

% li_IGA;
addpath(genpath('C:\Users\Administrator\Desktop\calfem-matlab-iga-master - (2)'));
%% ��ʼģ��(ֱ��)
main_IGA;  

%% �޸�ģ��
ANSWER1 = questdlg('Select the example to run.', ...
                  'chose the type of calculation', ...
                  'time_analysis', 'static_reanalysis_local ', 'static_reanalysis_global', 'defualt');
%% ����ʱ̬����            
if strcmpi(ANSWER1,  'time_analysis')    
    %%----����ʱ̬�ĳ���
    return;
end
%% �����ط����ľ�̬�ṹ�Ż�
if strcmpi(ANSWER1,   'static_reanalysis_local ')

%% �����ط������Ż� (local)
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
%% �����ط������Ż� (global)




