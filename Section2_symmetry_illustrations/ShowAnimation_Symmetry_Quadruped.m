function  ShowAnimation_Symmetry_Quadruped(T,Y_origin,P_origin,Y_mapped,P_mapped,color_plot,RecordKeyFrames)
% Plot settings
ScreenSize = get(0,'ScreenSize');

% Defining the size and positions of the plots
PlotSize      = [(2.0/10)*ScreenSize(3)   (12/16)*(2.0/10)*ScreenSize(3)];
PlotPositions = [(1/2)*ScreenSize(3)-(2/2 + 1/10)*PlotSize(1,1)  (2.0/10)*ScreenSize(4)+PlotSize(1,2)   PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(0/2 - 1/10)*PlotSize(1,1)  (2.0/10)*ScreenSize(4)+PlotSize(1,2)   PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(2/2 + 1/10)*PlotSize(1,1)  (1.5/10)*ScreenSize(4)                 PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(0/2 - 1/10)*PlotSize(1,1)  (1.5/10)*ScreenSize(4)                 PlotSize(1,1)  PlotSize(1,2)];

am = 'Detailed';
options = struct('AnimationMode',am);

%% Record Key Frames 
    if RecordKeyFrames == 1
%         tEvents = P(1:9);
%         tUpdate = sort(tEvents);
        tUpdate = sort([P_origin(1:9) T(round(length(T)/2))]);
        frameCount = 1;
        options.AnimationMode = 'Convinient';

        figure(105)
        graphOUTPUT_origin = SLIP_Animation_Quad(P_origin,PlotPositions(1,:),figure(105),options);
        for k = 1:length(tUpdate)
            y_origin = interp1(T' + linspace(0,1e-5,length(T)), Y_origin, tUpdate(k))';
            graphOUTPUT_origin.update(y_origin,P_origin,tUpdate(k)+0.001);
            fig = gcf;
            if  tUpdate(k) == T(round(length(T)/2))
                filen = 'f'+string(frameCount)+'_midstance_origin.pdf';
            elseif  tUpdate(k) == P_origin(9)
                filen = 'f'+string(frameCount)+'_apex_origin.pdf';
            else
                filen = 'f'+string(frameCount)+'_origin.pdf';
            end
            exportgraphics(fig,filen,'ContentType','vector')
            frameCount = frameCount + 1;
        end
        
        tUpdate = sort([P_mapped(1:9) T(round(length(T)/2))]);
        frameCount = 1;
        figure(205)
        graphOUTPUT_mapped = SLIP_Animation_Quad(P_mapped,PlotPositions(2,:),figure(205),options);
        for k = 1:length(tUpdate)
            y_mapped = interp1(T' + linspace(0,1e-5,length(T)), Y_mapped, tUpdate(k))';
            graphOUTPUT_mapped.update(y_mapped,P_mapped,tUpdate(k)+0.001);
            fig = gcf;
            if  tUpdate(k) == T(round(length(T)/2))
                filen = 'f'+string(frameCount)+'_midstance_mapped.pdf';
            elseif  tUpdate(k) == P_mapped(9)
                filen = 'f'+string(frameCount)+'_apex_mapped.pdf';
            else
                filen = 'f'+string(frameCount)+'_mapped.pdf';
            end
            exportgraphics(fig,filen,'ContentType','vector')
            frameCount = frameCount + 1;
        end
        
        return
        
    end
    
%% Plot Animation 

    % Plot Animation Before Mapping
    figure(105)
    graphOUTPUT_origin = SLIP_Animation_Quad(P_origin,PlotPositions(1,:),figure(105),options);
    graphOUTPUT_origin.Title.String = ["SLIP Model Animation","Before Mapping"];
    % Plot Periodic Orbit Before Mapping
    figure(106)
    poOUTPUT_origin    = SLIP_PeriodicOrbit_Quad(Y_origin,PlotPositions(2,:),figure(106),color_plot);
    poOUTPUT_origin.axes.Title.String = 'Periodic Orbit Before Mapping';


    % Plot Animation After Mapping
    figure(205)
    graphOUTPUT_mapped = SLIP_Animation_Quad(P_mapped,PlotPositions(3,:),figure(205),options);
    graphOUTPUT_mapped.Title.String = ["SLIP Model Animation","After Mapping"];
    % Plot Periodic Orbit After Mapping
    figure(206)
    poOUTPUT_mapped   = SLIP_PeriodicOrbit_Quad(Y_mapped,PlotPositions(4,:),figure(206),color_plot);
    poOUTPUT_mapped.axes.Title.String = 'Periodic Orbit After Mapping';

    n = round(T(end)*100); % # of frames per step
    tFrame = linspace(0, T(end), n+1);
    frameCount = 1;
    % If desired, every iteration a rendered picture is saved to disc.  This
    % can later be used to create a animation of the monopod.
    for j = 1:n

        y_origin    = interp1(T' + linspace(0,1e-5,length(T)), Y_origin,   tFrame(j));
        % Update Animation Before Mapping
        graphOUTPUT_origin.update(y_origin,P_origin,tFrame(j));
        % Update Periodic Orbit Before Mapping
        poOUTPUT_origin.update(y_origin);

        y_mapped    = interp1(T' + linspace(0,1e-5,length(T)), Y_mapped,   tFrame(j));
        % Update Animation After Mapping
        graphOUTPUT_mapped.update(y_mapped,P_mapped,tFrame(j));
        % Update Periodic Orbit After Mapping
        poOUTPUT_mapped.update(y_mapped);

        frameCount = frameCount + 1;
    end


end

