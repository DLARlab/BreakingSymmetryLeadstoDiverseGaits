function  ShowTrajectory_Quadruped_Demo(T,Y,GRF)

% Get screen setting: self tunning plot positions
ScreenSize = get(0,'ScreenSize');
PlotSize = [(2.5/10)*ScreenSize(3)   (12/16)*(2.5/10)*ScreenSize(3)
             ScreenSize(3)/5         (12/16)*ScreenSize(3)/5     ];
PlotPositions = [(1/2)*ScreenSize(3)-(3/2)*PlotSize(2,1)         (0.5/10)*ScreenSize(4)+PlotSize(2,2)                PlotSize(2,1)  PlotSize(2,2)
                 (1/2)*ScreenSize(3)-(1/2)*PlotSize(2,1)         (0.5/10)*ScreenSize(4)+PlotSize(2,2)                PlotSize(2,1)  PlotSize(2,2)
                 (1/2)*ScreenSize(3)+(1/2)*PlotSize(2,1)         (0.5/10)*ScreenSize(4)+PlotSize(2,2)                PlotSize(2,1)  PlotSize(2,2)];
%% Plot trajectories
    figure(101);set(gcf,'Position', PlotPositions(1,:));
    clf; grid on; plot(T,Y(:,2:6),'-');
    legend('dx','y','dy','phi','dphi');
    title('Trajectories of the Torso')
    
    figure(102);set(gcf,'Position', PlotPositions(2,:));
    clf; 
    subplot(2,1,1);
    plot(T,Y(:,7:10),'-');
    legend('alphaBL','dalphaBL','alphaFL','dalphaFL');
    grid on;
    title('Trajectories of the Left Legs')
    subplot(2,1,2);
    plot(T,Y(:,11:14),'-');
    legend('alphaBR','dalphaBR','alphaFR','dalphaFR');
    grid on;
    title('Trajectories of the Right Leg')

    [pks,locs] = findpeaks(Y(:,3));

%% Plot Ground Reaction Force
      
    FBLy = GRF(:,1);
    FFLy = GRF(:,2);
    FBRy = GRF(:,3);
    FFRy = GRF(:,4);
    
    
    grf = figure(203); clf; grid on; hold on;
    set(grf,'position',PlotPositions(3,:)); 
    plot(T,FBLy,'color',[217/256,83/256,25/256], 'LineWidth',2);  
    plot(T,FFLy,':','color',[217/256,83/256,25/256], 'LineWidth',2.5);  
    plot(T,FFRy,':','color',[0/256,114/256,189/256], 'LineWidth',2.5);
    plot(T,FBRy,'color',[0/256,114/256,189/256], 'LineWidth',2);
    xlim([0 T(end)]);
    title('Ground Reaction Forces')
%     ylim([0 3])

    legend('LH','LF','RF','RH'); 
    box on;
    xlabel('Stride Time  $[rad]$','Interpreter','LaTex','FontSize',12)
    ylabel('Vertical GRF  $[m_0g]$','Interpreter','LaTex','FontSize',12)
%     [pks,locs] = findpeaks(FFRy);
    M1 = max(FBLy);  M2 = max(FFLy);   M3 = max(FFRy);  M4 = max(FBRy);
    M5 = min(FBLy);  M6 = min(FFLy);   M7 = min(FFRy);  M8 = min(FBRy);
    if min([M5 M6 M7 M8])<0 
        ylim([min([M5 M6 M7 M8])*1.1  max([M1,M2,M3,M4])*1.1]);
    else
        ylim([0 max([M1,M2,M3,M4])*1.1]);
    end
    

end

