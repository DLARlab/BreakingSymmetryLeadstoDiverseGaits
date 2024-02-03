%% Section 2: exemplary solutions of pronking, bounding, half-bounding and galloping gaits
% The code in this section is used to show the solution multiple gaits.
% The code will show how Fig.3, Fig.5, Fig.6 are generated.
% Part 1: Searching for a single periodic solution
% Part 2: Search for a branch of periodic solution

% The readers can simply follow the instruction to finish this section.

close all
clear
clc


disp('...')
pause(2)
disp('Section 2: exemplary solutions illustration')
disp('...')

pause(2)
disp('The code in this section is used for animating different types of asymmetrical gaits mentioned in the paper.')
disp('...')

pause(2)
disp('And replicate Fig.3, Fig.5, Fig.6.')
disp('...')

pause(2)
disp('The readers can simply follow the instructions showing in the terminal.')
disp('...')


disp('(Hit Enter to Move on)')
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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol = scatter3(X(1),X(3),X(5),120,'filled','Parent',Roadmap.Children,...
                                           'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
pause(1)
current_sol_text =  text(X(1), X(3)-Roadmap.Children.YLim(2)*0.1, X(5),...
                     'Current Solution','HorizontalAlignment','center','Parent',Roadmap.Children,'FontSize',8);
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)+Roadmap.Children.YLim(2)*0.13, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.13, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.6*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.1, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.3*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.1, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.55*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.12, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show solution branch for half-bounding with hind leg spread and extended suspension
pause(2)
figure(1)
load('BD1_20_2_HE.mat')
[gait,abbr, color_plot, linetype] = Gait_Identification(results);
he = plot3(results(1,:),results(5,:),results(2,:),'LineWidth',3,'Color',color_plot,'LineStyle',linetype);
pause(1)
he_text = text(results(1,round(0.2*size(results,2))) - Roadmap.Children.XLim(2)*0.02, 0-Roadmap.Children.YLim(2)*0.35, 0, ...
    abbr,'Color',color_plot,'FontSize',12,'FontWeight','bold');
pause(1)

fprintf('\n') 
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,round(0.55*size(results,2)));
Para = results(23:end,round(0.3*size(results,2)));
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.12, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,96);
Para = results(23:end,96);
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)+Roadmap.Children.YLim(2)*0.2, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')

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
disp('Current gait of solution branch:')
pause(2)
disp(string(gait)+ ' gait')
pause(2)
disp("Abbreviate as '" + string(abbr) + "'")
pause(2)
disp('...')

% show current solution on the plot
X = results(1:22,115);
Para = results(23:end,115);
current_sol.XData = X(1);
current_sol.YData = X(5); 
current_sol.ZData = X(3); 
pause(1)
current_sol_text.Position =  [X(1), X(5)-Roadmap.Children.YLim(2)*0.1, X(3)];
pause(1)


pause(2)
disp('Showing the animation of exemplary solution...')
disp('...')

pause(2)
% simulation the system with initial condition
[residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para);
% show animation of the current solution
Animation = 1; PO = 0; RecordKeyFrames = 0;
ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)

pause(2)
disp('Animation Ends.')
disp('...')

disp('(Hit Enter to Move on)')
%pause()
fprintf('\n')
disp('...')


pause(2)
fprintf('\n') 
fprintf('\n')
disp('Section 2 complete.')
disp('...')
fprintf('\n') 


pause(2)
disp("Please open 'Section_3_Symmetry_Illustrations.m'.")