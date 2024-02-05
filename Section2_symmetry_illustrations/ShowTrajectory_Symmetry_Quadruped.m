function  ShowTrajectory_Symmetry_Quadruped(T,Y_origin,GRF_origin,Y_mapped,GRF_mapped)

% Get screen setting: self tunning plot positions
ScreenSize = get(0,'ScreenSize');
PlotSize = [(2.0/10)*ScreenSize(3)   (12/16)*(2.0/10)*ScreenSize(3)];
PlotPositions = [(1/2)*ScreenSize(3)-(3/2)*PlotSize(1,1)         (2.0/10)*ScreenSize(4)+PlotSize(1,2)                PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(1/2)*PlotSize(1,1)         (2.0/10)*ScreenSize(4)+PlotSize(1,2)                PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)+(1/2)*PlotSize(1,1)         (2.0/10)*ScreenSize(4)+PlotSize(1,2)                PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(3/2)*PlotSize(1,1)         (1.5/10)*ScreenSize(4)                              PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(1/2)*PlotSize(1,1)         (1.5/10)*ScreenSize(4)                              PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)+(1/2)*PlotSize(1,1)         (1.5/10)*ScreenSize(4)                              PlotSize(1,1)  PlotSize(1,2)];


    %% Plot trajectories after mapping
    figure(201);set(gcf,'Position', PlotPositions(4,:));
    clf; grid on; plot(T,Y_mapped(:,2:6),'-');
    legend('dx','y','dy','phi','dphi');
    title('Trajectories of the Torso After Mapping')
    
    figure(202);set(gcf,'Position', PlotPositions(5,:));
    clf; 
    subplot(2,1,1);
    plot(T,Y_mapped(:,7:10),'-');
    legend('alphaBL','dalphaBL','alphaFL','dalphaFL');
    grid on;
    title('Trajectories of the Left Legs After Mapping')
    subplot(2,1,2);
    plot(T,Y_mapped(:,11:14),'-');
    legend('alphaBR','dalphaBR','alphaFR','dalphaFR');
    grid on;
    title('Trajectories of the Right Leg After Mapping')

    [pks,locs] = findpeaks(Y_mapped(:,3));

%% Plot Ground Reaction Force after mapping
      
    FBLy = GRF_mapped(:,1);
    FFLy = GRF_mapped(:,2);
    FBRy = GRF_mapped(:,3);
    FFRy = GRF_mapped(:,4);
    
    
    GRF_mapped = figure(203); clf; grid on; hold on;
    set(GRF_mapped,'position',PlotPositions(6,:)); 
    plot(T,FBLy,'color',[217/256,83/256,25/256], 'LineWidth',2);  
    plot(T,FFLy,':','color',[217/256,83/256,25/256], 'LineWidth',2.5);  
    plot(T,FFRy,':','color',[0/256,114/256,189/256], 'LineWidth',2.5);
    plot(T,FBRy,'color',[0/256,114/256,189/256], 'LineWidth',2);
    xlim([0 T(end)]);
    title('Ground Reaction Forces After Mapping')
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
    
%% Plot trajectories before mapping
    figure(101);set(gcf,'Position', PlotPositions(1,:));
    clf; grid on; plot(T,Y_origin(:,2:6),'-');
    legend('dx','y','dy','phi','dphi');
    title('Trajectories of the Torso Before Mapping')
    
    figure(102);set(gcf,'Position', PlotPositions(2,:));
    clf; 
    subplot(2,1,1);
    plot(T,Y_origin(:,7:10),'-');
    legend('alphaBL','dalphaBL','alphaFL','dalphaFL');
    grid on;
    title('Trajectories of the Left Legs Before Mapping')
    subplot(2,1,2);
    plot(T,Y_origin(:,11:14),'-');
    legend('alphaBR','dalphaBR','alphaFR','dalphaFR');
    grid on;
    title('Trajectories of the Right Leg Before Mapping')

    [pks,locs] = findpeaks(Y_origin(:,3));

%% Plot Ground Reaction Force before mapping
      
    FBLy = GRF_origin(:,1);
    FFLy = GRF_origin(:,2);
    FBRy = GRF_origin(:,3);
    FFRy = GRF_origin(:,4);
    
    
    GRF_origin = figure(103); clf; grid on; hold on;
    set(GRF_origin,'position',PlotPositions(3,:)); 
    plot(T,FBLy,'color',[217/256,83/256,25/256], 'LineWidth',2);  
    plot(T,FFLy,':','color',[217/256,83/256,25/256], 'LineWidth',2.5);  
    plot(T,FFRy,':','color',[0/256,114/256,189/256], 'LineWidth',2.5);
    plot(T,FBRy,'color',[0/256,114/256,189/256], 'LineWidth',2);
    xlim([0 T(end)]);
    title('Ground Reaction Forces Before Mapping')
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

