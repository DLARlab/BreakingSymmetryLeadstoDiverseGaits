classdef SLIP_PeriodicOrbit_Quad < OutputCLASS 
        properties 
            fig; % The output window
            axes;
            % Patch objects used in the graphical representation
            Orbit;
            Poincare_Section;
            Current_Position;
            Text;
        end
      methods
        % Constructor:
        function obj = SLIP_PeriodicOrbit_Quad(Y,PlotPositions,fig,color_plot)
            obj.slowDown = 1;      % Run this in real time.
            obj.rate     = 0.05;   % with 25 fps
            
            obj.fig = fig;     clf(obj.fig);
            % Set window properties
            set(obj.fig, 'Name','Periodic Orbit');  % Window title
            set(obj.fig, 'Color','w');          % Background color
            set(obj.fig, 'Renderer','OpenGL');
            set(obj.fig, 'position', PlotPositions);
            
            obj.axes = axes(obj.fig);  hold on;
            obj.axes.View = [45 45];

            obj.axes.XLim = [min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1)   max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1)];
            obj.axes.YLim = [min(Y(:,4))*(1-sign(min(Y(:,4)))*0.1)   max(Y(:,4))*(1+sign(max(Y(:,4)))*0.1)];
            obj.axes.ZLim = [min(Y(:,6))-0.02  max(Y(:,6))+0.02];
            obj.axes.XLabel = xlabel('$\dot{q}_x  [\sqrt{gl_0}]$','Interpreter','LaTex','FontSize',15);  
            obj.axes.YLabel = ylabel('$\dot{q}_z  [\sqrt{gl_0}]$','Interpreter','LaTex','FontSize',15);
            obj.axes.ZLabel = zlabel('$\dot{q}_{pitch}  [rad/s]$','Interpreter','LaTex','FontSize',15);
            obj.axes.Title.String = 'Periodic Orbit of The Solution';


            obj.Orbit = plot3(Y(:,2),Y(:,4),Y(:,6),'LineWidth',2,'Color',[0 0 0],'Parent',obj.axes);
            
            % The edge of the section plane is self-adjustable: depending on min/max value, and their signs          
            obj.Poincare_Section = fill3([min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1) max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1) max(Y(:,2))*(1+sign(max(Y(:,2)))*0.1) min(Y(:,2))*(1-sign(min(Y(:,2)))*0.1)],[0 0 0 0],...
                                         [min(Y(:,6))-0.02 min(Y(:,6))-0.02 max(Y(:,6))+0.02 max(Y(:,6))+0.02],[0.7 0.7 0.7],...
                                         'Parent',obj.axes,'FaceAlpha',0.3);
            
            textPS = text( max(Y(:,2))*1.06, 0, (max(Y(:,6)) + min(Y(:,6)))*0.1,...
                           {'Poincare Section','$\dot{q}_z = 0$'},'Interpreter','LaTex','HorizontalAlignment','center',...
                           'Parent',obj.axes);
            textPO = text( min(Y(:,2))*0.98, min(Y(:,4)), (max(Y(:,6))+min(Y(:,6)))/2,...
                            'Periodic Orbit','HorizontalAlignment','center',...
                            'Parent',obj.axes);
            obj.Text = struct('textPS',textPS,'textPO',textPO);
            
            obj.Current_Position = scatter3(Y(1,2),Y(1,4),Y(1,6),120,'filled','Parent',obj.axes,...
                                           'MarkerEdgeColor',color_plot,'MarkerFaceColor',color_plot);
        end
        
        function obj = update(obj,y)
            figure(obj.fig)
            set(obj.Current_Position,'xData',y(2),'yData',y(4),'zData',y(6));
        end
      end
end