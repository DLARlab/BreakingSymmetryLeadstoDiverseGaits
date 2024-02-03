function  ShowAnimation_Quadruped(T,Y,P,color_plot,Animation,PO,RecordKeyFrames)
% Plot settings
ScreenSize = get(0,'ScreenSize');

% Defining the size and positions of the plots
PlotSize      = [(2.5/10)*ScreenSize(3)   (12/16)*(2.5/10)*ScreenSize(3)
                 ScreenSize(3)/5        (12/16)*ScreenSize(3)/5     ];
PlotPositions = [(1/2)*ScreenSize(3)-(2/2 + 1/10)*PlotSize(1,1)  (0.5/10)*ScreenSize(4)+PlotSize(2,2)  PlotSize(1,1)  PlotSize(1,2)
                 (1/2)*ScreenSize(3)-(0/2 - 1/10)*PlotSize(1,1)  (0.5/10)*ScreenSize(4)+PlotSize(2,2)  PlotSize(1,1)  PlotSize(1,2)];

am = 'Detailed';
options = struct('AnimationMode',am);

%% Plot Animation
    % Plot Animation
    if Animation == 1
        graphOUTPUT = SLIP_Animation_Quad(P,PlotPositions,options);
    end
    % Plot Periodic Orbit
    if PO == 1
        poOUTPUT    = SLIP_PeriodicOrbit_Quad(Y,PlotPositions,color_plot);
    end
    n = round(T(end)*100); % # of frames per step
    tFrame = linspace(0, T(end), n+1);
    frameCount = 1;
    % If desired, every iteration a rendered picture is saved to disc.  This
    % can later be used to create a animation of the monopod.
    for j = 1:n
        y    = interp1(T' + linspace(0,1e-5,length(T)), Y,   tFrame(j));
        if Animation == 1
            % Update Animation
            graphOUTPUT.update(y,P,tFrame(j));
        end
        if PO == 1
            % Update Periodic Orbit
            poOUTPUT.update(y);
        end
        frameCount = frameCount + 1;
    end

%% Record Key Frames
    if RecordKeyFrames == 1
%         tEvents = P(1:9);
%         tUpdate = sort(tEvents);
        tUpdate = sort([P(1:9) T(round(length(T)/2))]);
        frameCount = 1;
        options.AnimationMode = 'Convinient';
        graphOUTPUT = SLIP_Animation_Quad(P,PlotPositions,options);
        for k = 1:length(tUpdate)
            y = interp1(T' + linspace(0,1e-5,length(T)), Y, tUpdate(k))';
            graphOUTPUT.update(y,P,tUpdate(k)+0.001);
            fig = gcf;
            if  tUpdate(k) == T(round(length(T)/2))
                filen = 'f'+string(frameCount)+'_midstance.pdf';
            elseif  tUpdate(k) == P(9)
                filen = 'f'+string(frameCount)+'_apex.pdf';
            else
                filen = 'f'+string(frameCount)+'.pdf';
            end
            exportgraphics(fig,filen,'ContentType','vector')
            frameCount = frameCount + 1;
        end    
    end
end

