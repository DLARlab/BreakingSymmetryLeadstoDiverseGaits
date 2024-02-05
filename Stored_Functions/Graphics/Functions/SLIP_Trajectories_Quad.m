classdef SLIP_Trajectories_Quad < OutputCLASS 
    %SLIP_TRAJECTORIES_QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Torsofig;
        Torsoaxes;
        Torsoplot;
        Torsolegend;
        Legfig;
        Leftlegaxes;
        Leftlegplot;
        Leftleglegend;
        Rightlegtaxes;
        Rightlegtplot;
        Rightlegtlegend;
        
    end
    
    methods
        function obj = SLIP_Trajectories_Quad(T,Y,PlotPositions)
                    obj.slowDown = 1;      % Run this in real time.
                    obj.rate     = 0.004;   % with 250 fps
                    
                    % Trajectory Figure for Torso
                    obj.Torsofig = figure(105);
                    set(obj.Torsofig,'Position', PlotPositions(3,:)); 
                    clf(obj.Torsofig);% Set the plot size, position and clear it if there exist axes or plots on the current figure.
                    % Axes, trajectory lines and legend for Torso
                    obj.Torsoaxes = axes(obj.Torsofig); 
                    obj.Torsoaxes.XGrid = 'Off'; obj.Torsoaxes.YGrid = 'Off'; hold on; box(obj.Torsoaxes);
                    obj.Torsoaxes.XLabel = xlabel('Stride Time $[rad]$','Interpreter','LaTex','FontSize',12);
                    obj.Torsoaxes.YLabel = ylabel('Normalized States','Interpreter','LaTex','FontSize',12);
                    obj.Torsoaxes.Title.String = 'States and Velocities of the Torso';
                    obj.Torsoaxes.XLim = [0 T(end)];
                    obj.Torsoaxes.YLim = [min(min(Y(:,2:6)))*1.1  max(max(Y(:,2:6)))*1.1];

                    obj.Torsoplot = plot(T,Y(:,2:6),'-','Parent',obj.Torsoaxes);
                    obj.Torsolegend = legend('$\dot{x}$','$y$','$\dot{y}$','$\phi$','$\dot{\phi}$','Interpreter','LaTex');

                    % Trajectory Plot for Legs
                    obj.Legfig   = figure(106);
                    set(obj.Legfig,'Position', PlotPositions(4,:));
                    clf(obj.Legfig);
                    % Axes, trajectory lines and legend for Left legs
                    obj.Leftlegaxes = subplot(2,1,1);
                    obj.Leftlegaxes.XGrid = 'Off'; obj.Leftlegaxes.YGrid = 'Off'; hold on;
                    obj.Leftlegaxes.Title.String = 'States and Velocities of the Left Legs';
                    obj.Leftlegaxes.YLabel = ylabel('Normalized States','Interpreter','LaTex','FontSize',12);
                    obj.Leftlegaxes.XLim = [0 T(end)];
                    obj.Leftlegaxes.YLim = [min(min(Y(:,7:10)))*1.1  max(max(Y(:,7:10)))*1.1];

                    obj.Leftlegplot = plot(T,Y(:,7:10),'-');
                    obj.Leftleglegend = legend('$\alpha_{BL}$','$\dot{\alpha}_{BL}$','$\alpha_{FL}$','$\dot{\alpha}_{FL}$','Interpreter','LaTex');
                    
                
                    % Axes, trajectory lines and legend for right legs
                    obj.Rightlegtaxes = subplot(2,1,2);
                    obj.Rightlegtaxes.XGrid = 'Off'; obj.Rightlegtaxes.YGrid = 'Off'; hold on;
                    obj.Rightlegtaxes.Title.String = 'States and Velocities of the Right Legs';
                    obj.Rightlegtaxes.XLabel = xlabel('Stride Time $[rad]$','Interpreter','LaTex','FontSize',12);
                    obj.Rightlegtaxes.YLabel = ylabel('Normalized States','Interpreter','LaTex','FontSize',12);
                    obj.Rightlegtaxes.XLim = [0 T(end)];
                    obj.Rightlegtaxes.YLim = [min(min(Y(:,11:14)))*1.1  max(max(Y(:,11:14)))*1.1];

                    obj.Rightlegtplot = plot(T,Y(:,11:14),'-');
                    obj.Rightlegtlegend = legend('$\alpha_{BR}$','$\dot{\alpha}_{BR}$','$\alpha_{FR}$','$\dot{\alpha}_{FR}$','Interpreter','LaTex');
        end
        
        function obj = update(obj,T_,Y_)
            
                    set(obj.Torsoplot(1),   'XData',T_ ,'YData',Y_(:,2))
                    set(obj.Torsoplot(2),   'XData',T_ ,'YData',Y_(:,3))
                    set(obj.Torsoplot(3),   'XData',T_ ,'YData',Y_(:,4))
                    set(obj.Torsoplot(4),   'XData',T_ ,'YData',Y_(:,5))
                    set(obj.Torsoplot(5),   'XData',T_ ,'YData',Y_(:,6))
                  
                    set(obj.Leftlegplot(1), 'XData',T_ ,'YData',Y_(:,7))
                    set(obj.Leftlegplot(2), 'XData',T_ ,'YData',Y_(:,8))
                    set(obj.Leftlegplot(3), 'XData',T_ ,'YData',Y_(:,9))                  
                    set(obj.Leftlegplot(4), 'XData',T_ ,'YData',Y_(:,10))
                    
                    set(obj.Rightlegtplot(1),'XData',T_, 'YData',Y_(:,11))
                    set(obj.Rightlegtplot(2),'XData',T_, 'YData',Y_(:,12))
                    set(obj.Rightlegtplot(3),'XData',T_, 'YData',Y_(:,13))
                    set(obj.Rightlegtplot(4),'XData',T_, 'YData',Y_(:,14))
        end
    end
end

