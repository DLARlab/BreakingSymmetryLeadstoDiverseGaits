%% Section 2: Illustration of Temporal Symmetry, Time-Reversal Symmetry, and Morphological Symmetry
% This code illustrates the symmetries mentioned in the paper.
% The code consists of three parts:
% Part 1: Temporal symmetry
% Part 2: Time-reversal symmetry
% Part 3: Morphological symmetry

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

disp('...')
disp('Section 2: illustration of symmetries')
disp('...')

pause(2)
disp('The code in this section is used for illustrating different types of symmetries mentioned in the paper.')
disp('...')


pause(2)
disp('The readers can simply follow the instructions showing in the terminal.')
disp('...')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: illustration of temporal symmetry
pause(2)
disp('Part 1: Temporal symmetry')
disp('...')
pause(2)
disp('The solution is symmetric under the evolution after a certain time T.')
disp('...')
pause(2)
disp('This means the solution will return to itself after simulating the system for time T.')
disp('...')
pause(2)
disp('This property is retained no matter where the evolution starts along the periodic orbit.')
disp('...')
fprintf('\n')


pause(2)
load('BD1_20_2_FG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
disp('Loading solution...')
disp('...')

pause(2)
% load a solution
X_origin = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
% simulation the system with initial condition
[residual_origin,T_origin,Y_origin,P_origin,GRF_origin,Y_EVENT_origin] = Quadrupedal_ZeroFun_v2(X_origin,Para);
residual_origin(end) = []; % remove the residual for Poincare section
disp('Simulating the system with the initial condition...')
disp('...')
pause(2)
disp('X_origin(t+T) - X_origin(t) = ' + string(norm(residual_origin)))
disp('...')


% pick another point on the periodic orbit
pause(2)
X_permutation = []';
X_permutation(1:13)  = Y_origin(round(0.3*length(T_origin)),2:14);
X_permutation(14:21) = X_origin(14:21) - T_origin(round(0.3*length(T_origin)));
X_permutation(22)    = X_origin(22);  
X_permutation        = X_permutation';
disp('Picking another point along the orbit...')
disp('...')

% simulation the system from another point on the periodic orbit
pause(2)
[residual_mapped,T_mapped,Y_mapped,P_mapped,GRF_mapped,Y_EVENT_mapped] = Quadrupedal_ZeroFun_v2(X_permutation,Para);
residual_mapped(end) = []; % remove the residual for Poincare section
disp('Simulating the system from another point along the orbit...')
disp('...')
pause(2)
disp('X_temporal(t+T) - X_temporal(t) = ' + string(norm(residual_mapped)))
disp('...')
% make sure the two simulations are in the same length
if length(T_mapped) < length(T_origin)
    Y_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), Y_origin,   T_mapped);
    GRF_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), GRF_origin,   T_mapped);
    T = T_mapped;
elseif length(T_mapped) > length(T_origin)
    Y_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), Y_mapped,   T_origin);
    GRF_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), GRF_mapped,   T_origin);
    T = T_origin;
else
    T = T_origin;
end


pause(2)
ShowTrajectory_Symmetry_Quadruped(T,Y_origin,GRF_origin,Y_mapped,GRF_mapped)
disp('Plotting the trajectories of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')

pause(2)
RecordKeyFrames = 0;
ShowAnimation_Symmetry_Quadruped(T,Y_origin,P_origin,Y_mapped,P_mapped,color_plot,RecordKeyFrames)
disp('Plotting the animations of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: illustration of time-reversal symmetry
clear

disp('Part 2: Time-reversal symmetry')
disp('...')
pause(2)
disp("The solution remains temporally symmetric when the system's velocities are reversed and the configuration is unchanged.")
disp('...')
pause(2)
disp('Meaning the solution will also return to itself after simulate the system for the same time period T.')
disp('...')
pause(2)
disp('This property is retained no matter where the evolution starts along the periodic orbit.')
disp('...')
fprintf('\n')


pause(2)
load('BD1_20_2_BG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
disp('Loading solution...')
disp('...')

pause(2)
% load a solution
X_origin = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
% simulation the system with initial condition
[residual_origin,T_origin,Y_origin,P_origin,GRF_origin,Y_EVENT_origin] = Quadrupedal_ZeroFun_v2(X_origin,Para);
residual_origin(end) = []; % remove the residual for Poincare section
disp('Simulating the system with the initial condition...')
disp('...')
pause(2)
disp('X_origin(t+T) - X_origin(t) = ' + string(norm(residual_origin)))
disp('...')



pause(2)
X_permutation = X_origin;
% reverse all the velocities but leave the configuration unchanges
X_permutation(1:2:13)  = -X_origin(1:2:13);
% flip the touchdown/liftoff event and shift the timings
X_permutation([15 17 19 21 14 16 18 20]) = X_permutation(22) - X_permutation([14 16 18 20 15 17 19 21]);
disp('Reversing the velocities of the solution...')
disp('...')
pause(2)
disp('               [ -1  0  0  0  0  0  0  0  0  0  0  0  0  ]  [ d(q_x)/dt     ]')
disp('               [  0  1  0  0  0  0  0  0  0  0  0  0  0  ]  [ q_y           ]')
disp('               [  0  0 -1  0  0  0  0  0  0  0  0  0  0  ]  [ d(q_y)/dt     ]')
disp('               [  0  0  0  1  0  0  0  0  0  0  0  0  0  ]  [ q_pitch       ]')
disp('               [  0  0  0  0 -1  0  0  0  0  0  0  0  0  ]  [ d(q_pitch)/dt ]')
disp('               [  0  0  0  0  0  1  0  0  0  0  0  0  0  ]  [ q_BL          ]')
disp('  \varphi(x) = [  0  0  0  0  0  0 -1  0  0  0  0  0  0  ]  [ d(q_BL)/dt    ]')
disp('               [  0  0  0  0  0  0  0  1  0  0  0  0  0  ]  [ q_FL          ]')
disp('               [  0  0  0  0  0  0  0  0 -1  0  0  0  0  ]  [ d(q_FL)/dt    ]')
disp('               [  0  0  0  0  0  0  0  0  0  1  0  0  0  ]  [ q_BR          ]')
disp('               [  0  0  0  0  0  0  0  0  0  0 -1  0  0  ]  [ d(q_BR)/dt    ]')
disp('               [  0  0  0  0  0  0  0  0  0  0  0  1  0  ]  [ q_FR          ]')
disp('               [  0  0  0  0  0  0  0  0  0  0  0  0 -1  ]  [ d(q_FR)/dt    ]')
disp('...')
pause(5)

% simulation the system with reversed velocities and the same configuration
pause(2)
[residual_mapped,T_mapped,Y_mapped,P_mapped,GRF_mapped,Y_EVENT_mapped] = Quadrupedal_ZeroFun_v2(X_permutation,Para);
residual_mapped(end) = []; % remove the residual for Poincare section
disp('Simulating the system with reversed velocities and the same configuration...')
disp('...')
pause(2)
disp('X_time_reversal(t+T) - X_time_reversal(t) = ' + string(norm(residual_mapped)))
disp('...')
% make sure the two simulations are in the same length
if length(T_mapped) < length(T_origin)
    Y_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), Y_origin,   T_mapped);
    GRF_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), GRF_origin,   T_mapped);
    T = T_mapped;
elseif length(T_mapped) > length(T_origin)
    Y_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), Y_mapped,   T_origin);
    GRF_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), GRF_mapped,   T_origin);
    T = T_origin;
else
    T = T_origin;
end

pause(2)
ShowTrajectory_Symmetry_Quadruped(T,Y_origin,GRF_origin,Y_mapped,GRF_mapped)
disp('Plotting the trajectories of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')

pause(2)
RecordKeyFrames = 0;
ShowAnimation_Symmetry_Quadruped(T,Y_origin,P_origin,Y_mapped,P_mapped,color_plot,RecordKeyFrames)
disp('Plotting the animations of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')

pause(2)
disp('The reader can also compare the key frames to validate the time-reversal symmetry.')
disp('...')
pause(2)
disp('The time sequence of both the oringinal solution and the time-reversed solution will be recorded.')
disp('...')
pause(2)
disp("To do so, please enter 'Y' in the terminal, or hit any button to move on.")
disp('...')

prompt = 'Recode key frames? Yes[Y]/No(press any button).';
input = input(prompt,"s");
if string(input)=='Y'
    RecordKeyFrames = 1;
    ShowAnimation_Symmetry_Quadruped(T,Y_origin,P_origin,Y_mapped,P_mapped,color_plot,RecordKeyFrames)
end
fprintf('\n')
disp('...')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: illustration of morphological symmetry
clear

disp('Part 3: Morphological symmetry')
disp('...')
pause(2)
disp('The solution is symmetric under the permutation of the states of the legs.')
disp('...')
pause(2)
disp('The trajectories will be identical to the trajectories before permutation.')
disp('...')
pause(2)
disp('This property is retained no matter where the evolution starts along the periodic orbit.')
disp('...')
fprintf('\n')


pause(2)
load('BD1_20_2_HE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
disp('Loading solution...')
disp('...')

pause(2)
% load a solution
X_origin = results(1:22,round(0.55*size(results,2)));
Para = results(23:end,round(0.55*size(results,2)));
% simulation the system with initial condition
[residual_origin,T_origin,Y_origin,P_origin,GRF_origin,Y_EVENT_origin] = Quadrupedal_ZeroFun_v2(X_origin,Para);
residual_origin(end) = []; % remove the residual for Poincare section
disp('Simulating the system with the initial condition...')
disp('...')
pause(2)
disp('X_origin(t+T) - X_origin(t) = ' + string(norm(residual_origin)))
disp('...')



pause(2)
X_permutation = X_origin;
% permutate the states of the legs
X_permutation([8 9 12 13 16 17 20 21])  = X_origin([12 13 8 9 20 21 16 17]);
disp('Permutating the states of the legs...')
disp('...')
pause(2)
disp('                      [  1  0  0  0  0  0  0  0  0  0  0  0  0  ]  [ d(q_x)/dt     ]')
disp('                      [  0  1  0  0  0  0  0  0  0  0  0  0  0  ]  [ q_y           ]')
disp('                      [  0  0  1  0  0  0  0  0  0  0  0  0  0  ]  [ d(q_y)/dt     ]')
disp('                      [  0  0  0  1  0  0  0  0  0  0  0  0  0  ]  [ q_pitch       ]')
disp('                      [  0  0  0  0  1  0  0  0  0  0  0  0  0  ]  [ d(q_pitch)/dt ]')
disp('                      [  0  0  0  0  0  1  0  0  0  0  0  0  0  ]  [ q_BL          ]')
disp('  \sigma_{FL,FR}(x) = [  0  0  0  0  0  0  1  0  0  0  0  0  0  ]  [ d(q_BL)/dt    ]')
disp('                      [  0  0  0  0  0  0  0  0  0  0  0  1  0  ]  [ q_FL          ]')
disp('                      [  0  0  0  0  0  0  0  0  0  0  0  0  1  ]  [ d(q_FL)/dt    ]')
disp('                      [  0  0  0  0  0  0  0  0  0  1  0  0  0  ]  [ q_BR          ]')
disp('                      [  0  0  0  0  0  0  0  0  0  0  1  0  0  ]  [ d(q_BR)/dt    ]')
disp('                      [  0  0  0  0  0  0  0  1  0  0  0  0  0  ]  [ q_FR          ]')
disp('                      [  0  0  0  0  0  0  0  0  1  0  0  0  0  ]  [ d(q_FR)/dt    ]')
disp('...')
pause(5)

% simulation the system with permutated leg states
pause(2)
[residual_mapped,T_mapped,Y_mapped,P_mapped,GRF_mapped,Y_EVENT_mapped] = Quadrupedal_ZeroFun_v2(X_permutation,Para);
residual_mapped(end) = []; % remove the residual for Poincare section
disp('Simulating the system with permutated leg states...')
disp('...')
pause(2)
disp('X_permutation(t+T) - X_permutation(t) = ' + string(norm(residual_mapped)))
disp('...')
% make sure the two simulations are in the same length
if length(T_mapped) < length(T_origin)
    Y_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), Y_origin,   T_mapped);
    GRF_origin = interp1(T_origin' + linspace(0,1e-5,length(T_origin)), GRF_origin,   T_mapped);
    T = T_mapped;
elseif length(T_mapped) > length(T_origin)
    Y_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), Y_mapped,   T_origin);
    GRF_mapped = interp1(T_mapped' + linspace(0,1e-5,length(T_mapped)), GRF_mapped,   T_origin);
    T = T_origin;
else
    T = T_origin;
end

pause(2)
ShowTrajectory_Symmetry_Quadruped(T,Y_origin,GRF_origin,Y_mapped,GRF_mapped)
disp('Plotting the trajectories of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')

pause(2)
RecordKeyFrames = 0;
ShowAnimation_Symmetry_Quadruped(T,Y_origin,P_origin,Y_mapped,P_mapped,color_plot,RecordKeyFrames)
disp('Plotting the animations of the solution before and after mapping.')
disp('...')
disp('(Hit Any Button to Move on)')
pause()
close all
fprintf('\n')
disp('...')


pause(2)
fprintf('\n') 
fprintf('\n')
disp('Section 2 complete.')
disp('...')
fprintf('\n') 
