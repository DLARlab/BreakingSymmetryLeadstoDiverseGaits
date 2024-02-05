%% Section 3: Example of Solution Searching Process
% The code in this section is used to show the solution searching process of a simplified model.
% The section consists of two parts:
% Part 1: Searching for a single periodic solution
% Part 2: Search for a branch of periodic solution

% Readers can simply follow the instructions to complete this section.

close all
clear
clc

restoredefaultpath
current_path = matlab.desktop.editor.getActiveFilename;
current_path = current_path(1:strfind(current_path,'BreakingSymmetryLeadstoDiverseGaits')+35); 

disp('...')
disp('Adding project path...')
disp('...')
addpath(genpath(string(current_path)+'Section1_exemplary_solutions'))
addpath(genpath(string(current_path)+'Section2_symmetry_illustrations'))
addpath(genpath(string(current_path)+'Section3_solution_searching'))
addpath(genpath(string(current_path)+'Stored_Functions'))
disp('Path Initialized.')
disp('...')

fprintf('\n') 
disp('...')
pause(2)
disp('Section 3: example of solution searching process')
disp('...')

pause(2)
disp('This code demonstrates the solution searching process for a simplified model.')
disp('...')

pause(2)
disp('Readers can simply follow the instructions shown in the terminal.')
disp('...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part1: Searching for a single periodic solution
pause(2)
fprintf('\n') 
disp('Execution begins. Follow the steps as outlined below:')

pause(2)
disp('...')
disp('Part 1: Searching for a single periodic solution')

pause(3)
disp('...')
disp("Current Stage:  Giving a random guess for the system's states... ")
disp('...')


% Giving an random guess of the system states
pause(2)
X_initial = [1.3          % d(q_x)/dt
             1.1          % q_y
             0            % d(q_y)/dt
             0            % q_pitch
             0            % d(q_pitch)/dt
             0            % q_BL
             1.0          % d(q_BL)/dt
             0            % q_FL
             1.0          % d(q_FL)/dt
             0            % q_BR
             1.0          % d(q_BR)/dt
             0            % q_FR
             1.0 ];       % d(q_FR)/dt
         
disp('X_initial = [1.3          % d(q_x)/dt')
disp('             1.1          % q_y')
disp('             0            % d(q_y)/dt')
disp('             0            % q_pitch')
disp('             0            % d(q_pitch)/dt')
disp('             0            % q_BL')
disp('             1.0          % d(q_BL)/dt')
disp('             0            % q_FL')
disp('             1.0          % d(q_FL)/dt')
disp('             0            % q_BR')
disp('             1.0          % d(q_BR)/dt')
disp('             0            % q_FR')
disp('             1.0          % d(q_FR)/dt')
disp('...')

% Giving an initial guess of the eventing timings
pause(2) 
Timing_initial = [0.6         % t_BL_TD
                  0.9         % t_BL_LO
                  0.6         % t_FL_TD    
                  0.9         % t_FL_LO
                  0.6         % t_BR_TD
                  0.9         % t_BR_LO
                  0.6         % t_FR_TD
                  0.9         % t_FR_LO
                  1.3 ];      % t_Apex(Stride Period)
fprintf('\n') 
disp('Giving a random guess of the event timings for the solution:') 
fprintf('\n')                              
disp('Timing_initial = [0.6         % t_BL_TD')
disp('                  0.9         % t_BL_LO')
disp('                  0.6         % t_FL_TD  ')
disp('                  0.9         % t_FL_LO')
disp('                  0.6         % t_BR_TD')
disp('                  0.9         % t_BR_LO')
disp('                  0.6         % t_FR_TD')
disp('                  0.9         % t_FR_LO')
disp('                  1.3 ];      % t_Apex(Stride Period)')
disp('...')

pause(2)
solution_initial = [X_initial;Timing_initial];
disp('solution_initial = [X_initial;Timing_initial]')
disp('...') 

pause(2)
disp('A detailed definition of the solution can be found in Stored_Functions\System_Dynamics\Quadrupedal_ZeroFun_v2.m.')
disp('...')

disp('(Hit Any Button to Move on)')
pause()
fprintf('\n') 
disp('...')

% define the parameters of the system
Para = [10                  % k: linear leg stiffness of legs
        20                  % ks: swing stiffness of legs, omega = sqrt(ks);
        2                   % J: torso pitching inertia
        1                   % l: ratio between resting leg length and main body l/L, usually set to be 1.
        0                   % osa: Resting angle of swing leg motion 
        0.5                 % lb: distance from COM to hip joint/length of torso
        1];                 % kr: ratio of linear stiffness between back and front legs: kr = kb/kf
    
% Description of how the dynamic function helps to find solutions with symmetries.
pause(2)
disp('Code Description:')
fprintf('\n') 
disp('The dynamic function Quadrupedal_ZeroFun_v2.m integrates the system with given initial conditions and timing variables.')
pause(2)
fprintf('\n') 
disp("While evoluting, the function will also measure the 'symmetry error', which is supposed to be zero if the symmetry is retained. ")
pause(2)
fprintf('\n') 
disp("The solution search is equivalent to solving the zero problem for 'symmetry error'. ")
disp('...')
fprintf('\n') 

pause(2)
disp('Example:')
pause(1)
fprintf('\n')
disp('Solving for the solution with temporal symmetry (periodic solution).')
fprintf('\n') 
disp("'symmetry error' = x_end - x_initial.")
pause(3)
fprintf('\n') 
disp('The zero problem will be solved using fsolve.m from MATLAB.' )
disp('...')
fprintf('\n')


pause(2)
disp('Starting solving the zero problem.')
pause(2)
% options for solving zero function     
numOPTS = optimset('Algorithm','levenberg-marquardt',... 
                   'ScaleProblem','jacobian',...
                   'Display','iter',...
                   'MaxFunEvals',15000,...
                   'MaxIter',3000,...
                   'UseParallel', false,...
                   'TolFun',1e-9,...
                   'TolX',1e-12);
% solve the zero problem            
[X1_real_solution, ~] = fsolve(@(X) Quadrupedal_ZeroFun_v2(X,Para), solution_initial, numOPTS);
pause(2)
fprintf('\n') 
disp('Zero problem solved.')
disp('...')
pause(2)

% identifying the type of the solution(gait identification)
[gait,abbr, color_plot, linetype] = Gait_Identification(X1_real_solution);
fprintf('\n') 
disp('Gait identified as:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
pause(2)
disp('...')


pause(3)
disp('Showing the animation of current solution.')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X1_real_solution,Para);
% show animation of the current solution
Animation = 1; PO = 1; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)



disp('(Hit Any Button to Move on)')
pause()
close all
disp('...')
disp('Showing the trajectories of current solution.')
disp('...')
% show trajectories of the current solution
ShowTrajectory_Quadruped_Demo(T,Y,GRF)


pause(5)
disp('Part 1 complete.')
disp('...')
fprintf('\n') 

disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part2: Searching for a branch of periodic solutions
pause(2)
fprintf('\n') 
disp('...')
disp('Part 2: Searching for a branch of periodic solutions')

pause(2)
disp('...')
disp('Code Description:')
fprintf('\n') 
disp('This part shows the solution searching process for a entire branch using NumericalContinuation1D_Quadruped_Demo.m. ')
pause(2)
fprintf('\n') 
disp('The algorithm implements a prediction-corrector method to search for the solution branch.' )
pause(2)
fprintf('\n') 
disp('To start with, the algorithm requires two solutions X1 X2 as inputs, and will predict the new solution using the difference between the two solutions:')
pause(2)
fprintf('\n') 
disp('X_predict = X2 + (X2 - X1)')
pause(2)
fprintf('\n') 
disp("The algorithm will then solve for a zero probelm mentioned in Part 1 to 'correct' the predicted solution.")
pause(2)
fprintf('\n') 
disp("The entire solution branch can be found by iterations.")
disp('...')
fprintf('\n') 


disp('(Hit Any Button to Move on)')
pause()
fprintf('\n') 
disp('...')

% give a random prediction for the new initial guess
pause(2)
X2_predict = X1_real_solution + 0.05*(rand-0.5);
disp('Providing an initial guess for the second solution.')
pause(2)
disp('...')
disp('X2_predict = X_real_solution + 0.05*(rand-0.5)')
disp('...')
fprintf('\n')

pause(5)
disp('Solve for the zero problem')
disp('...')
[X2_real_solution, ~] = fsolve(@(X) Quadrupedal_ZeroFun_v2(X,Para), X2_predict, numOPTS);
pause(2)
fprintf('\n') 
disp('Zero problem solved.')
disp('...')


pause(2)
fprintf('\n') 
disp('Running the numerical continuation algorithm to search for the entire solution branch')
pause(2)
disp('...')
disp('The demo algorithm will terminate after finding 10 solutions.')
pause(2)
disp('...')
disp('To search for the entire branch, please change the number of iterations in NumericalContinuation1D_Quadruped_Demo.m')
disp('...')
pause(2)
% options for the algorithm
numOPTS = optimset('Algorithm','levenberg-marquardt',...
                   'ScaleProblem','jacobian',...
                   'Display','iter',...
                   'MaxFunEvals',50000,...
                   'MaxIter',3000,...
                   'UseParallel', false,...
                   'TolFun',1e-9,...
                   'TolX',1e-12);
radius = 0.05;
[results,flag] = NumericalContinuation1D_Quadruped_Demo(X1_real_solution,X2_real_solution,Para,radius,numOPTS);

pause(2)
disp('...')
disp('Demo algorithm terminated.')
disp('...')
fprintf('\n') 

pause(2)
disp('Part 2 complete.')
disp('...')
fprintf('\n') 

pause(2)
fprintf('\n') 
fprintf('\n')
disp('Section 3 complete.')
disp('...')
fprintf('\n') 




pause(2)
disp("Demo code complete! Thank you for following along!")
