% *************************************************************************
% classdef SLIP_Model_Graphics(p) < OutputCLASS
%
% Two dimensional graphics of a SLIP model.
%
% The graphics object must be initialized with the vector of system
% parameters p.
%
%
% Properties: - NONE
% Methods:    - NONE
%
%
% Created by C. David Remy on 07/10/2011
% MATLAB 2010a - Windows - 64 bit
%
% Documentation:
%  'A MATLAB Framework For Gait Creation', 2011, C. David Remy (1), Keith
%  Buffinton (2), and Roland Siegwart (1),  International Conference on
%  Intelligent Robots and Systems, September 25-30, San Francisco, USA 
%
% (1) Autonomous Systems Lab, Institute of Robotics and Intelligent Systems, 
%     Swiss Federal Institute of Technology (ETHZ) 
%     Tannenstr. 3 / CLA-E-32.1
%     8092 Zurich, Switzerland  
%     cremy@ethz.ch; rsiegwart@ethz.ch
%
% (2) Department of Mechanical Engineering, 
%     Bucknell University
%     701 Moore Avenue
%     Lewisburg, PA-17837, USA
%     buffintk@bucknell.edu
%
%   See also OUTPUTCLASS.
%
classdef SLIP_Animation_Quad< OutputCLASS 
    % Private attributes:
    properties 
        fig;
        axes;
        Body;  
        Leg_BL;
        Leg_FL;
        Leg_BR;
        Leg_FR;
        COM;
        Ground;
        PhaseDiagram;
        Title;
        options
    end
    % Public methods:
    methods
        % Constructor:
        function obj = SLIP_Animation_Quad(P,PlotPositions,fig,options)
            obj.slowDown = 1;      % Run this in real time.
            obj.rate     = 0.05;   % with 25 fps
            obj.options  = options;
            
            % Copy the parameter vector: COM location
            lb     = P(15);
            l_leg  = P(16);
            
            % Initialize the graphics and set properties
            obj.fig = fig; clf(obj.fig);
            set(obj.fig, 'Name','SLIP model');  % Window title
            set(obj.fig, 'Color','w');          % Background color
            set(obj.fig, 'Renderer','OpenGL');
            set(obj.fig, 'Position', PlotPositions);

            % Initialize the axes and set properties: axis off
            obj.axes = axes(obj.fig);
            outerpos = obj.axes.OuterPosition;
            ti = obj.axes.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            set(obj.axes, 'Position',[left bottom ax_width ax_height],'Box','Off')
%             obj.axes.Position = [left bottom ax_width ax_height];
            obj.axes.XAxis.Visible = 'Off'; 
            obj.axes.YAxis.Visible = 'Off';


            
            % Define some arbitrary states:
            x        = 1;
            y        = 1.2;
            phi_body = 0;
            
            % Right Legs
            obj.Leg_BL = DrawLegs([x-lb;y],l_leg, 0.3);
            obj.Leg_FL = DrawLegs([x+(1-lb);y],l_leg,-0.1); 
            
            % Draw Ground
            obj.Ground = DrawGround; 
         
            % Main Body
            obj.Body  = DrawBody([x,y,phi_body],lb);                      
            %Centor of Mass
            if ~(lb == 0.5)
                radius = 0.15;
                obj.COM = DrawCOM([x,y,phi_body],radius);
            end
            
            % Left Legs
            obj.Leg_BR = DrawLegs([x-lb;y],l_leg, 0.3);
            obj.Leg_FR = DrawLegs([x+(1-lb);y],l_leg,-0.1);
            
            % PhaseDigram
            if string(obj.options.AnimationMode) == 'Detailed'
                obj.Title = text(x-1.3,1.9,'SLIP Model Animation','FontSize',15,'FontWeight','Bold');
                obj.PhaseDiagram = DrawPhaseDiagram(x,P);
            end
           

        end   
        % Updated function.Is called by the integrator:
        function obj = update(obj,y,P,T)
            
            figure(obj.fig)
            [LegLength, LegAngle, BodyJPos, BackJPos, FrontJPos] = ComputeJoint_LegLA(y,P,T);
                     
            % Left Legs
            SetLegs(BackJPos,LegLength.BL,LegAngle.BL,obj.Leg_BL);
            SetLegs(FrontJPos,LegLength.FL,LegAngle.FL,obj.Leg_FL);
            % Main
            SetBody(BodyJPos, obj.Body,P(15));
            % Centor of Mass
            if ~(P(15) == 0.5)
                radius = 0.15;
                SetCOM(BodyJPos,obj.COM,radius);
            end         
            % Right Legs
            SetLegs(BackJPos,LegLength.BR,LegAngle.BR,obj.Leg_BR);
            SetLegs(FrontJPos,LegLength.FR,LegAngle.FR,obj.Leg_FR);        
            % phase diagram
            
            if string(obj.options.AnimationMode)=='Detailed'
                set(obj.Title,'Position',[BodyJPos(1)-1.3 1.9]);
                SetPhaseDiagram(BodyJPos(1),P,obj.PhaseDiagram)
            end
            
            axis([-1.5 + BodyJPos(1),1.5 + BodyJPos(1),-0.1,2]);
%             axis([-2 ,10,-0.2,2]);
            drawnow();
        end
    end
end

%% Functions that draw the parts of the body

% Create patches of main body in the Constructor. Save handles for update function.
function  BodyH = DrawBody(BodyJPos,lb)

    [ x1, y1, f,v ]  = ComputeBodyGraphics(BodyJPos,lb);

    % White background (face color white, no line)
    b1  = patch('XData',x1,'YData',y1,'LineStyle','none','FaceColor',[1 1 1]); 
    % Shade lines (Grey lines)
    b2  = patch('faces', f, 'vertices', v,...
           'linewidth',3,'FaceColor',[1 1 1],'EdgeColor',0.8*[1 1 1]); 
    % Black Outline (No face color white, black outline)
    b3  = patch('XData',x1,'YData',y1,'linewidth',4,'FaceColor','none');

    % Save all the handles for update function
    BodyH = struct('B_bg',b1,'B_sha',b2,'B_out',b3);

end


% Create patches of legs in the Constructor. Save handles for update function.
function  LegParts=DrawLegs(vecS,l_leg, gamma_leg)

    [LegVertices, LegFaces] = ComputeLegGraphics(vecS,l_leg,gamma_leg);
    % Spring Part 1************************************************************
    L1   = patch('faces', LegFaces.L_Sp1, 'vertices', LegVertices.L_Sp1,...
           'linewidth',5,'EdgeColor',[245 131 58]/256); % color blue

    % Upper Leg****************************************************************
    % 1. outline of upper leg
    L2_1  = patch('faces',LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO,...
          'LineStyle','none','FaceColor',[1 1 1]);
    % 2. Draw shaded region in the upper leg
    L3   = patch('faces', LegFaces.L_Ups,  'vertices', LegVertices.L_Ups,...
          'linewidth',3,'FaceColor','none','EdgeColor',0.8*[1 1 1]); % grey
    % 3. Outline of upper leg  
    L2_2  = patch('faces',LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO,...
          'linewidth',3,'FaceColor','none');

    % Lower Leg ***************************************************************
    L4  = patch('faces', LegFaces.L_low, 'vertices', LegVertices.L_low,...
          'linewidth',3,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);

    % Spring Part 2 ***********************************************************
    L5   = patch('faces',LegFaces.L_Sp2, 'vertices', LegVertices.L_Sp1,...
          'linewidth',5,'EdgeColor',[245 131 58]/256);

    % Save all the handles for update function
    LegParts = struct('L_Sp1',L1,'L_UpB',L2_1,'L_Upo',L2_2,'L_Ups',L3,'L_low',L4,'L_Sp2',L5);

end

% Create patches of Centor of mass in the Constructor. Save handles for update function.
function COMH = DrawCOM(BodyJPos,radius)

    RotM = [ cos(BodyJPos(3)), -sin(BodyJPos(3));
             sin(BodyJPos(3)),  cos(BodyJPos(3))];
    % Draw Center of Mass
    alpha = linspace(0, pi*2, 40);
    % vert_x_out = sin(alpha)*0.2;
    % vert_y_out = cos(alpha)*0.2;
    vert_out = [sin(alpha)*radius*0.5
                cos(alpha)*radius*0.5];
    vert_out = RotM*vert_out + BodyJPos(1:2)'*ones(1,size(vert_out,2));


    b1 = patch(vert_out(1,:), vert_out(2,:),'white','linewidth',5); 

    alpha = linspace(0, pi/2, 10);
    vert = [0,sin(alpha)*radius*0.75,0
                 0,cos(alpha)*radius*0.75,0];
    vert_x = [vert(1,:);vert(1,:);-vert(1,:);-vert(1,:)];
    vert_y = [vert(2,:);-vert(2,:);-vert(2,:);vert(2,:)];
    Vert = zeros(8,size(vert,2));
    for i = 1:4
        Vert(2*i-1:2*i,:) = RotM*[vert_x(i,:);vert_y(i,:)];
    end

    vert_x = [Vert(1,:);Vert(3,:);Vert(5,:);Vert(7,:)]' + BodyJPos(1);
    vert_y = [Vert(2,:);Vert(4,:);Vert(6,:);Vert(8,:)]' + BodyJPos(2);
    b2 = patch(vert_x, vert_y, cat(3,[1 0 1 0], [1 0 1 0],[1 0 1 0]),'LineWidth',3);


    COMH = struct('B_COMOuter',b1,'B_COMInner',b2);

end

% Create patches of Ground in the Constructor. Save handles for update function.
function GroundH = DrawGround

        % Draw the ground. It reaches from -2.5 to +6.5.
        h   = 0.01; % Height of the bar at the top
        n   = 5000;  % Number of diagonal stripes in the shaded area
        s   = 0.05; % Spacing of the stripes
        w   = 0.01; % Width of the stripes
        ext = 0.1;  % Length of the stripes

        % Create vertices by shifting a predefined pattern 'n' times to the right:
        v = [     -50,0;
            repmat([     0,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([  -ext,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([-ext+w,-ext-h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            repmat([     w,    -h],n,1) + [-50+linspace(0,s*n,n)',zeros(n,1)];
            -50+s*n+w,0];
        % Connect to faces:
        f = [1,2,4*n+1,4*n+2;
             repmat([0,n,2*n,3*n],n,1) + repmat((1:n)',1,4)+1];
         
        vert_x_out = [-15 100 100 -15];
        vert_y_out = [0 0 -20 -20];
        ground = patch(vert_x_out, vert_y_out,'white');   
        lines  = patch('faces', f, 'vertices', v);
        
        GroundH = struct('ground',ground,'lines',lines);
end

% Create patches of Phase Diagram in the Constructor. Save handles for update function.
function PhaseDiagram = DrawPhaseDiagram(x,P)
    %Compute the position of vertices
    [vertices_box,pos_text,pos_PhaseBars] = ComputePhaseDiagram(x,P);

    % A box will be drawn to include the phase bars and discriptions
    % Define vertices of box 
    box_x = vertices_box(1,:);
    box_y = vertices_box(2,:);
    phase_box = patch(box_x, box_y,'white','LineWidth',1.5);


    % Show the descriptions for the phase bars
    t_bl = text(pos_text(1,1),pos_text(1,2),'LH','FontSize',12.5);
    t_fl = text(pos_text(2,1),pos_text(2,2),'LF','FontSize',12.5);
    t_fr = text(pos_text(3,1),pos_text(3,2),'RF','FontSize',12.5);
    t_br = text(pos_text(4,1),pos_text(4,2),'RH','FontSize',12.5);

    % Drawing the phase bars: draw two layers, one for stance face and one for flight phase.
    % If tTD<tLO, draw flight phase(white) and then cover it with stance phase(black);
    % If tTD>tLO, draw stance phase(white) and then cover it with flight phase(black).

    % Define the color of phase bar, currently black
    cp = [0 0 0];

    % Event timings that determine the sequence and the length of phase bars
    pos_BL = pos_PhaseBars(1:4,:);
    pos_BR = pos_PhaseBars(5:8,:);
    pos_FL = pos_PhaseBars(9:12,:);
    pos_FR = pos_PhaseBars(13:16,:);
    tEvents = P(1:9);
    % Back left phase bar
    if tEvents(1)<tEvents(2)
        phasebl_ = patch(pos_BL(1,:), pos_BL(2,:), 'white','EdgeColor','none');
        phasebl  = patch(pos_BL(3,:), pos_BL(4,:), cp,'EdgeColor','none');
    elseif tEvents(1)>tEvents(2)
        phasebl_ = patch(pos_BL(1,:), pos_BL(2,:), cp,'EdgeColor','none');
        phasebl  = patch(pos_BL(3,:), pos_BL(4,:), 'white','EdgeColor','none');
    end

    % Back right phase bar
    if tEvents(5)<tEvents(6)
        phasebr_ = patch(pos_BR(1,:), pos_BR(2,:),'white','EdgeColor','none');
        phasebr  = patch(pos_BR(3,:), pos_BR(4,:), cp,'EdgeColor','none');
    elseif tEvents(5)>tEvents(6)
        phasebr_ = patch(pos_BR(1,:), pos_BR(2,:), cp,'EdgeColor','none');
        phasebr  = patch(pos_BR(3,:), pos_BR(4,:), 'white','EdgeColor','none');
    end

    % Front left phase bar
    if tEvents(3)<tEvents(4)
        phasefl_ = patch(pos_FL(1,:), pos_FL(2,:),'white','EdgeColor','none');
        phasefl  = patch(pos_FL(3,:), pos_FL(4,:), cp,'EdgeColor','none');
    elseif tEvents(3)>tEvents(4)
        phasefl_ = patch(pos_FL(1,:), pos_FL(2,:), cp,'EdgeColor','none');
        phasefl  = patch(pos_FL(3,:), pos_FL(4,:), 'white','EdgeColor','none');
    end


    % Front right phase bar
    if tEvents(7)<tEvents(8)
        phasefr_ = patch(pos_FR(1,:), pos_FR(2,:),'white','EdgeColor','none');
        phasefr  = patch(pos_FR(3,:), pos_FR(4,:), cp,'EdgeColor','none');
    elseif tEvents(7)>tEvents(8)
        phasefr_ = patch(pos_FR(1,:), pos_FR(2,:), cp,'EdgeColor','none');
        phasefr  = patch(pos_FR(3,:), pos_FR(4,:), 'white','EdgeColor','none');
    end


    PhaseDiagram = struct('Phase_bl',phasebl,'Phase_br',phasebr,'Phase_fl',phasefl,'Phase_fr',phasefr,...
                          'Phase_bl_',phasebl_,'Phase_br_',phasebr_,'Phase_fl_',phasefl_,'Phase_fr_',phasefr_,...
                          'Phase_box',phase_box,'Text_bl',t_bl,'Text_fl',t_fl,'Text_fr',t_fr,'Text_br',t_br);

end

%% Set Functions in updating the animation
% Set the patches of main body when figure is updated.
function SetBody(BodyJPos, Body,lb)

    [ x1, y1, f,v ] = ComputeBodyGraphics(BodyJPos,lb);

    set(Body.B_bg , 'xData', x1, 'yData',    y1);                   
    set(Body.B_sha, 'faces',  f, 'vertices',  v);
    set(Body.B_out, 'xData', x1, 'yData',    y1);   

end

% Set the patches of legs when figure is updated.
function  SetLegs(vecS,l_leg, gamma_leg,Leghandle)

    [LegVertices, LegFaces] = ComputeLegGraphics(vecS,l_leg,gamma_leg);
    % Spring Part 1 (Zigzag line)**********************************************
    set(Leghandle.L_Sp1,'faces', LegFaces.L_Sp1,  'vertices', LegVertices.L_Sp1);
    % Upper Leg****************************************************************
    % 1. Background of upper leg
    set(Leghandle.L_UpB,'faces', LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO );
    % 2. Shade the upper leg
    set(Leghandle.L_Ups,'faces', LegFaces.L_Ups,  'vertices', LegVertices.L_Ups);
    % 3. Outline the upper leg
    set(Leghandle.L_Upo,'faces', LegFaces.L_UpBO, 'vertices', LegVertices.L_UpBO );
    % Lower Leg ***************************************************************
    set(Leghandle.L_low,'faces', LegFaces.L_low,  'vertices', LegVertices.L_low );
    % Spring Part 2 (parallel)*************************************************
    set(Leghandle.L_Sp2,'faces', LegFaces.L_Sp2,  'vertices', LegVertices.L_Sp1);

end

% Set the patches of center of mass when figure is updated.
function SetCOM(BodyJPos,COM,radius)
    RotM = [ cos(BodyJPos(3)), -sin(BodyJPos(3));
             sin(BodyJPos(3)),  cos(BodyJPos(3))];
    % Draw Center of Mass
    alpha = linspace(0, pi*2, 40);
    % vert_x_out = sin(alpha)*0.2;
    % vert_y_out = cos(alpha)*0.2;
    vert_out = [sin(alpha)*radius*0.5
                cos(alpha)*radius*0.5];
    vert_out = RotM*vert_out + BodyJPos(1:2)'*ones(1,size(vert_out,2));

    alpha = linspace(0, pi/2, 10);
    vert = [0,sin(alpha)*radius*0.75,0
                 0,cos(alpha)*radius*0.75,0];
    vert_x = [vert(1,:);vert(1,:);-vert(1,:);-vert(1,:)];
    vert_y = [vert(2,:);-vert(2,:);-vert(2,:);vert(2,:)];
    Vert = zeros(8,size(vert,2));
    for i = 1:4
        Vert(2*i-1:2*i,:) = RotM*[vert_x(i,:);vert_y(i,:)];
    end

    vert_x = [Vert(1,:);Vert(3,:);Vert(5,:);Vert(7,:)]' + BodyJPos(1);
    vert_y = [Vert(2,:);Vert(4,:);Vert(6,:);Vert(8,:)]' + BodyJPos(2);

    set(COM.B_COMOuter, 'xData', vert_out(1,:), 'yData', vert_out(2,:))
    set(COM.B_COMInner, 'xData', vert_x, 'yData', vert_y)
end


% Set the patches of phase diagram when figure is updated.
function SetPhaseDiagram(x,P,PhaseDiagram)
    %Compute the position of vertices
    [vertices_box,pos_text,pos_PhaseBars] = ComputePhaseDiagram(x,P);

    % Define vertices of box 
    box_x = vertices_box(1,:);
    box_y = vertices_box(2,:);
    set(PhaseDiagram.Phase_box,'xData',box_x,'yData', box_y)


    % Define the position of text
    set(PhaseDiagram.Text_bl,'Position',pos_text(1,:))
    set(PhaseDiagram.Text_fl,'Position',pos_text(2,:))
    set(PhaseDiagram.Text_fr,'Position',pos_text(3,:))
    set(PhaseDiagram.Text_br,'Position',pos_text(4,:))

    % Event timings that determine the sequence and the length of phase bars
    pos_BL = pos_PhaseBars(1:4,:);
    pos_BR = pos_PhaseBars(5:8,:);
    pos_FL = pos_PhaseBars(9:12,:);
    pos_FR = pos_PhaseBars(13:16,:);

    % Back left phase bar
    set(PhaseDiagram.Phase_bl_,'xData',pos_BL(1,:),'yData', pos_BL(2,:))
    set(PhaseDiagram.Phase_bl,'xData', pos_BL(3,:),'yData', pos_BL(4,:))

    % Back right phase bar
    set(PhaseDiagram.Phase_br_,'xData',pos_BR(1,:),'yData', pos_BR(2,:))
    set(PhaseDiagram.Phase_br,'xData', pos_BR(3,:),'yData', pos_BR(4,:))

    % Front left phase bar
    set(PhaseDiagram.Phase_fl_,'xData',pos_FL(1,:),'yData', pos_FL(2,:))
    set(PhaseDiagram.Phase_fl,'xData', pos_FL(3,:),'yData', pos_FL(4,:))


    % Front right phase bar
    set(PhaseDiagram.Phase_fr_,'xData',pos_FR(1,:),'yData', pos_FR(2,:))
    set(PhaseDiagram.Phase_fr,'xData', pos_FR(3,:),'yData', pos_FR(4,:))

end

