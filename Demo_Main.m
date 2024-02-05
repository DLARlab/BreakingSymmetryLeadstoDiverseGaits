%% This file illustrates the exemplary solutions of a quadrupdel system: pronking, bounding, half-bounding and galloping.
% The symmetries of the solutions are also demonstrated and illustrated, as well as the solution search process.
% The code consistes of three sections:
% Section 1: Exemplary Solutions
% Section 2: Symmetry Illustrations
% Section 3: Solution Searching

% Readers can simply hit the 'Run' button to follow the whole demo code.
% Readers can also run any individual section if interested by hitting the 'Run Section' button.


%% Section 1: Exemplary Solutions for Pronking, Bounding, Half-Bounding, and Galloping Gaits
% This code illustrates solutions for multiple gaits.
% The code will show how Fig.3, Fig.5, Fig.6 are generated.

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
disp('Section 1: Illustration of Exemplary Solutions')
disp('...')

pause(2)
disp('This code is for animating the different types of asymmetrical gaits mentioned in the paper.')
disp('...')

pause(2)
disp('And replicate Fig.3, Fig.5, Fig.6.')
disp('...')

pause(2)
disp('Readers can simply follow the instructions to complete this section.')
disp('...')


%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')



% initialize figure
Roadmap = figure(1);
% plot settings
ScreenSize = get(0,'ScreenSize');
PlotPositionSize = [(1/2)*ScreenSize(3)-(0/2 - 1/10)*(2.5/10)*ScreenSize(3)  (0.5/10)*ScreenSize(4)+(12/16)*ScreenSize(3)/5 ...
                    (2.5/10)*ScreenSize(3)                                    (12/16)*(2.5/10)*ScreenSize(3)];

set(gcf,'Position', PlotPositionSize)
box on; hold on
xlabel('$\dot{x}  [\sqrt{gl_0}]$','Interpreter','LaTex','FontSize',15)
ylabel('$\dot{\varphi}  [rad]$','Interpreter','LaTex','FontSize',15)
zlabel('$y [l_0]$','Interpreter','LaTex','FontSize',15)

xlim([0 10])          ;     xticks(0:2:10)          ;
ylim([-0.05 0.05])    ;     yticks(-0.05:0.025:0.05);

title('Fig.3  Pronking / Bounding Solutions')
view([0 90])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for pronking 
pause(2)
load('PK_20_2.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
pk = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
pk_text = text(results(1,round(0.05*size(results,2))), results(3,round(0.05*size(results,2))) + Roadmap.Children.YLim(2)*0.1, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol = scatter3(X(1),X(3),X(5),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(3)-Roadmap.Children.YLim(2)*0.1, X(5),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(2)


% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('The animation ends.')
disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for bounding with gathered suspension 
pause(2)
figure(1)
load('BD1_20_2_BG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
bg = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
bg_text = text(results(1,round(0.2*size(results,2))) - Roadmap.Children.XLim(2)*0.1, results(5,round(0.2*size(results,2))) + Roadmap.Children.YLim(2)*0.1, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)+Roadmap.Children.YLim(2)*0.13, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)


% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for bounding with Extended suspension 
pause(2)
figure(1)
load('BD1_20_2_BE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
be = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
be_text = text(results(1,round(0.2*size(results,2))) - Roadmap.Children.XLim(2)*0.13, results(5,round(0.2*size(results,2))) - Roadmap.Children.YLim(2)*0.08, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.13, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)


% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for half-bounding with front leg spread and gathered suspension
pause(2)
figure(1)
xlim([2 7]);  xticks(2:1:7)  
title('Fig.5  Half-Bounding Solutions')

transparency = 0.2;
pk.Color = 1 - transparency * (1 - pk.Color);
pk_text.Color = 1 - transparency * (1 - pk_text.Color);
pk_text.Position = pk_text.Position + [2 0 0];

bg.Color = 1 - transparency * (1 - bg.Color);
bg_text.Color = 1 - transparency * (1 - bg_text.Color);
bg_text.Position = bg_text.Position + [2.5  0.01 0];
be.Color = 1 - transparency * (1 - be.Color);
be_text.Color = 1 - transparency * (1 - be_text.Color);
be_text.Position = be_text.Position + [2.5 -0.01 0];

load('BD1_20_2_FG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
fg = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
fg_text = text(results(1,round(0.3*size(results,2))) - Roadmap.Children.XLim(2)*0.03, results(5,round(0.3*size(results,2))) + Roadmap.Children.YLim(2)*0.20, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.6*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.1, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)



disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for half-bounding with hind leg spread and gathered suspension
pause(2)
figure(1)
load('BD1_20_2_HG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
hg = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
hg_text = text(results(1,round(0.2*size(results,2))) - Roadmap.Children.XLim(2)*0.05, results(5,round(0.2*size(results,2))) + Roadmap.Children.YLim(2)*0.15, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.13, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)




disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for half-bounding with front leg spread and extended suspension
pause(2)
figure(1)
load('BD1_20_2_FE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
fe = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
fe_text = text(results(1,round(0.2*size(results,2))) - Roadmap.Children.XLim(2)*0.10, 0-Roadmap.Children.YLim(2)*0.25, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.55*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.12, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)




disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for half-bounding with hind leg spread and extended suspension
pause(2)
figure(1)
load('BD1_20_2_HE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
he = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
he_text = text(results(1,round(0.2*size(results,2))) + Roadmap.Children.XLim(2)*0.02, 0-Roadmap.Children.YLim(2)*0.35, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.55*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.13, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)



% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for galloping with gathered suspension
pause(2)
figure(1)

fg.Color = 1 - transparency * (1 - fg.Color);
fg_text.Color = 1 - transparency * (1 - fg_text.Color);
hg.Color = 1 - transparency * (1 - hg.Color);
hg_text.Color = 1 - transparency * (1 - hg_text.Color);

fe.Color = 1 - transparency * (1 - fe.Color);
fe_text.Color = 1 - transparency * (1 - fe_text.Color);
he.Color = 1 - transparency * (1 - he.Color);
he_text.Color = 1 - transparency * (1 - he_text.Color);

title('Fig.6  Galloping Solutions')


load('BD1_20_2_GG.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
gg = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
gg_text = text(results(1,round(0.5*size(results,2))) + Roadmap.Children.XLim(2)*0.05, 0+Roadmap.Children.YLim(2)*0.35, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,96);
Para = results(23:end,96);

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)+Roadmap.Children.YLim(2)*0.20, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)



% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')

% Reset the position for the point(current solution)
delete(current_sol)
delete(current_sol_text)
close(figure(205))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for galloping with gathered suspension
pause(2)
figure(1)

load('BD1_20_2_GE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
ge = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
ge_text = text(results(1,round(0.5*size(results,2))) + Roadmap.Children.XLim(2)*0.10, 0-Roadmap.Children.YLim(2)*0.35, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current Gait of the Solution Branch:')
% pause(2)
disp(string(gait)+ ' gait')
% pause(2)
disp("Abbreviated as '" + string(abbr) + "'")
% pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,115);
Para = results(23:end,115);

current_sol = scatter3(X(1),X(5),X(3),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(5)-Roadmap.Children.YLim(2)*0.1, X(3),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)



% pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')
disp('Please do not change the position of the animation window.')
disp('...')


% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped_Demo(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
% disp('The animation ends.')
% disp('...')

%disp('(Hit Any Button to Move on)')
%pause()
fprintf('\n')
disp('...')


pause(2)
fprintf('\n') 
fprintf('\n')
disp('Section 1 complete.')
disp('...')
fprintf('\n') 


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
