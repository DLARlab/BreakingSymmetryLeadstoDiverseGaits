%% Section 0: initialization
% Initialization of the code, please hit the 'Run' button.
close all
clear
clc

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
addpath(genpath(current_path+slash+'Section1_exemplary_solutions'))
addpath(genpath(current_path+slash+'Section2_symmetry_illustrations'))
addpath(genpath(current_path+slash+'Section3_solution_searching'))

addpath(genpath(current_path+slash+'Stored_Functions')) 
pause(2)
disp('Section 0 complete')
disp('...')

pause(2)
disp("Please open 'Section_1_Exemplary_Solutions.m' and simply hit the 'Run' button.")