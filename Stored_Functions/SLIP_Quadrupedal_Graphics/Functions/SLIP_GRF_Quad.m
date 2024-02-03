classdef SLIP_GRF_Quad < OutputCLASS 
    %SLIP_GRF_QUAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties       
        fig;
        axes;
        BLPlot
        FLPlot
        BRPlot
        FRPlot
        legend;
    end
    
    methods
        function obj = SLIP_GRF_Quad(T,GRF,PlotPositions)
            
                    obj.slowDown = 1;      % Run this in real time.
                    obj.rate     = 0.004;   % with 250 fps
                    
                    FBL = GRF(:,1);
                    FFL = GRF(:,2);
                    FBR = GRF(:,3);
                    FFR = GRF(:,4);
                    % Find the peak values to determine the limit for axis
                    M1 = max(FBL);  M2 = max(FFL);   M3 = max(FFR);  M4 = max(FBR);
                    M5 = min(FBL);  M6 = min(FFL);   M7 = min(FFR);  M8 = min(FBR);
                    
                    % Figure for GRF plotting
                    obj.fig = figure(107); clf(obj.fig);
                    set(obj.fig,'position',PlotPositions(5,:)); 
                    % Axes for GRF plotting
                    obj.axes = axes(obj.fig);
                    obj.axes.XGrid = 'On'; obj.axes.YGrid = 'On'; hold on; box(obj.axes);
                    obj.axes.XLim = [0 T(end)];
                    if min([M5 M6 M7 M8])<0 
                        obj.axes.YLim = [min([M5 M6 M7 M8])*1.1  max([M1,M2,M3,M4])*1.1];
                    else
                        obj.axes.YLim = [0 max([M1,M2,M3,M4])*1.1];
                    end
                    obj.axes.XLabel = xlabel('Stride Time  $[rad]$','Interpreter','LaTex','FontSize',12);
                    obj.axes.YLabel = ylabel('Vertical GRF  $[m_0g]$','Interpreter','LaTex','FontSize',12);
                    obj.axes.Title.String = 'Ground Reaction Forces Applying on Legs'; 

                    obj.BLPlot = plot(T,FBL,'color',[217/256,83/256,25/256], 'LineWidth',2);  
                    obj.FLPlot = plot(T,FFL,':','color',[217/256,83/256,25/256], 'LineWidth',2.5);  
                    obj.FRPlot = plot(T,FFR,':','color',[0/256,114/256,189/256], 'LineWidth',2.5);
                    obj.BRPlot = plot(T,FBR,'color',[0/256,114/256,189/256], 'LineWidth',2);

                    obj.legend = legend('LH','LF','RF','RH'); 

        end
        
        function obj = update(obj,T_,GRF_)
                    FBL_ = GRF_(:,1);
                    FFL_ = GRF_(:,2);
                    FBR_ = GRF_(:,3);
                    FFR_ = GRF_(:,4);
                    set(obj.BLPlot,'XData',T_,'YData',FBL_)
                    set(obj.FLPlot,'XData',T_,'YData',FFL_)
                    set(obj.BRPlot,'XData',T_,'YData',FBR_)
                    set(obj.FRPlot,'XData',T_,'YData',FFR_)
        end
    end
end

