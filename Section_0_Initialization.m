%% Section 0: initialization
close all
clear
clc


fprintf('\n') 
disp('...')
pause(2)
disp('Section 0 : Enviroment Initializtion')
disp('...')

pause(2)
restoredefaultpath
current_path = string(pwd);
% slash and back slash
if ispc
    slash = '\';
else
    slash = '/';
end
disp('Adding project path...')
disp('...')
addpath(genpath(current_path+slash+'Section1_solution_searching'))
addpath(genpath(current_path+slash+'Section2_exemplary_solutions'))
addpath(genpath(current_path+slash+'Section3_symmetry_illustrations'))

addpath(genpath(current_path+slash+'Stored_Functions')) 
pause(2)
disp('Section complete')
disp('...')

pause(2)
disp("Please open 'Section_1_Solution_Searching.m'.")