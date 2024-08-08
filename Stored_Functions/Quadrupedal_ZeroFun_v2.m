%% This version of the dynamic model is designed for tunable system parameters: Para = [k,ks,J,lb,l]
%  The code is using ode function to integrate the system dynamics.
%  The script cnsists of:  1. Defining initial states, parameters, and phase
%                          2. Integrate the system using different ODE function base the phase (timing variables)
%                          3. The ODE functions: Stance Dynamics and Swing Dynamics
%                          4. Functions that compute the stance leg accelerations
%                          5. Computing the ground reaction forces
%  In this version v2, the stance dynamics are calculated based on input lb.

%  The output of the script consists of: 1. Residual values defined by periodic and holonomic constraints
%                                        2. States, Time, Parameters of the system
%                                        3. Vertical Ground Reaction Force
%                                        4. The States at Touchdown and Liftoff Timings

%% Dynamics Start from here: 
function [residual,T,Y,P,GRFs,Y_EVENT] = Quadrupedal_ZeroFun_v2(X,Para)
 
    %**********************************************************************
    % Parameter Preparation
    %**********************************************************************
    % Define model parameters:
    M  = 1;       % Set the torso mass to be 1
    
    % Identical Parameters
    k   = Para(1); % linear leg stiffness of legs
    ks  = Para(2); % swing stiffness of legs, omega = sqrt(ks);
    J   = Para(3); % torso pitching inertia
    l   = Para(4); % Ratio between resting leg length and main body l/L, usually set to be 1.
    osa = Para(5); % Resting angle of swing leg motion 
    
    % Asymmetrical Parameters
    lb = Para(6); % distance from COM to hip joint/length of torso
    kr = Para(7); % ratio of linear stiffness between back and front legs: kr = kb/kf, when kb+kf = 2*k;

    kb = 2*k/(1+ 1/kr); % linear stifness of back legs
    kf = 2*k/(1+ kr);   % linear stifness of front legs
    
    
    %**********************************************************************
    % States Preparation: the model has 7 DOFs, 14 states in total. The
    % initial condition also contains 8 timing variables, used to predefine
    % the touchdown/liftoff event. The dynamics equations are selected
    % based on the events defined.
    %**********************************************************************
    % Define the initial states of the integration
    N = 14; % 10 continuous states (including x)
    x0        = 0;
    dx0       = X(1); 
    y0        = X(2);
    dy0       = X(3);
    phi0      = X(4);
    dphi0     = X(5);    
    alphaBL0  = X(6);   % Left legs
    dalphaBL0 = X(7);
    alphaFL0  = X(8);
    dalphaFL0 = X(9);
    
    alphaBR0  = X(10);  % Right legs
    dalphaBR0 = X(11);
    alphaFR0  = X(12);
    dalphaFR0 = X(13);    
    %%%%%%%%%%%%%%%%%
    % Event timing:
    tBL_TD    = X(14);
    tBL_LO    = X(15);
    tFL_TD    = X(16);
    tFL_LO    = X(17);
    
    tBR_TD    = X(18);
    tBR_LO    = X(19);
    tFR_TD    = X(20);
    tFR_LO    = X(21);
    
    tAPEX     = X(22);



    %**********************************************************************
    % Event Timing Regulations: ensure 0 < t_i < t_APEX
    %**********************************************************************
    tBL_TD_ = tBL_TD;
    if  tBL_TD < 0
         tBL_TD_ = tBL_TD + tAPEX;
    end
    if tBL_TD > tAPEX
         tBL_TD_ = tBL_TD - tAPEX;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tBR_TD_ = tBR_TD;
    if tBR_TD < 0
         tBR_TD_ = tBR_TD + tAPEX;
    end
    if tBR_TD > tAPEX
         tBR_TD_ = tBR_TD - tAPEX;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tBL_LO_ = tBL_LO;
    if tBL_LO < 0
         tBL_LO_ = tBL_LO + tAPEX;
    end
    if tBL_LO > tAPEX
         tBL_LO_ = tBL_LO - tAPEX;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tBR_LO_ = tBR_LO;
    if tBR_LO < 0
         tBR_LO_ = tBR_LO + tAPEX;
    end
    if tBR_LO > tAPEX
         tBR_LO_ = tBR_LO - tAPEX;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tFL_LO_ = tFL_LO;
    if tFL_LO < 0
         tFL_LO_ = tFL_LO + tAPEX;
    end
    if tFL_LO > tAPEX
         tFL_LO_ = tFL_LO - tAPEX;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while tFR_LO < 0
         tFR_LO = tFR_LO + tAPEX;
    end
    while tFR_LO > tAPEX
         tFR_LO = tFR_LO - tAPEX;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tFL_TD_ = tFL_TD;
    if tFL_TD < 0
         tFL_TD_ = tFL_TD + tAPEX;
    end
    if tFL_TD > tAPEX
         tFL_TD_ = tFL_TD - tAPEX;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tFR_TD_ = tFR_TD;
    if tFR_TD < 0
         tFR_TD_ = tFR_TD + tAPEX;
    end
    if tFR_TD > tAPEX
         tFR_TD_ = tFR_TD - tAPEX;
    end




    %**********************************************************************
    % Integration Scheme: Chop the integration process into 8 pieces by using the event timing; the dynamic equations are determined by the timing variables
    %**********************************************************************
    % Set up start of integration:
    T_START = 0;
    Y_START = [x0, dx0, y0, dy0, phi0, dphi0,...
               alphaBL0, dalphaBL0, alphaFL0, dalphaFL0,...
               alphaBR0, dalphaBR0, alphaFR0, dalphaFR0];
    % Integrate motion in 8 steps, which are determined by the order of the event times (iEVENT(i) is the Eventnumber of the ith event)
    [tEVENT,iEVENT] = sort([tBL_TD_,tBL_LO_,tFL_TD_,tFL_LO_,tBR_TD_,tBR_LO_,tFR_TD_,tFR_LO,tAPEX]);
    % Prepare output:
    T = [];
    Y = [];
    Y_EVENT = zeros(9,N);
    
    for i = 1:9 %Integrate motion i/5
        % Figure out the current contact configuration (this is used in the
        % dynamics function)
        t_ = (T_START+tEVENT(i))/2;
        

        % Event detection conditions: determine the dynamic equations used for integration
        if ((t_>tBL_TD_ && t_<tBL_LO_ && tBL_TD_<tBL_LO_) || ((t_<tBL_LO_ || t_>tBL_TD_) && tBL_TD_>tBL_LO_))
            contactBL = true;
        else
            contactBL = false;
        end
        if ((t_>tFL_TD_ && t_<tFL_LO_ && tFL_TD_<tFL_LO_) || ((t_<tFL_LO_ || t_>tFL_TD_) && tFL_TD_>tFL_LO_))
            contactFL = true;
        else
            contactFL = false;
        end
        
        
        if ((t_>tBR_TD_ && t_<tBR_LO_ && tBR_TD_<tBR_LO_) || ((t_<tBR_LO_ || t_>tBR_TD_) && tBR_TD_>tBR_LO_))
            contactBR = true;
        else
            contactBR = false;
        end
        if ((t_>tFR_TD_ && t_<tFR_LO && tFR_TD_<tFR_LO) || ((t_<tFR_LO || t_>tFR_TD_) && tFR_TD_>tFR_LO))
            contactFR = true;
        else
            contactFR = false;
        end
        
        % Set up solver     
        %************************
        % Variable time step solver:
        % Setup ode solver options: 
            
%         options = odeset('RelTol',1e-12,'AbsTol',1e-12);
        
            if  abs(T_START - tEVENT(i))<1e-12 % T_START == tEVENT(i)  
                % make sure the time interval is valid
                Y_PART = Y_START;
                T_PART = T_START;
            else
                if X(1) < 15
                    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
                    [T_PART,Y_PART] = ode45(@ode,[T_START,tEVENT(i)],Y_START,options);
                else % change solver conditions if the Froude number is larger than 15
                    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
                    [T_PART,Y_PART] = ode45(@ode,[T_START,tEVENT(i)],Y_START,options);
                end
            end   
        
        % Event handlers:
        if iEVENT(i)==1
            % If this is EVENT 1, append left hind touchdown BL:
            T_PART=[T_PART;T_PART(end)]; 
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv1 = [Y_PART(end,[1 3 5 7 9]), Y_PART(end,[1 3 5 7 9]+1), zeros(1,5), lb]';
            Y_PART(end,8) = Func_alphaB_VA_v2(Pv1); % velocity reset at touchdown
        end
        
        if iEVENT(i)==3
            % If this is EVENT 3, append left front touchdown FL:
            T_PART=[T_PART;T_PART(end)];
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv2 = [Y_PART(end,[1 3 5 7 9]), Y_PART(end,[1 3 5 7 9]+1), zeros(1,5), lb]';
            Y_PART(end,10) = Func_alphaF_VA_v2(Pv2); % velocity reset at touchdown
        end
        
        if iEVENT(i)==5
            % If this is EVENT 5, append right hind touchdown BR:
            T_PART=[T_PART;T_PART(end)]; 
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv3 = [Y_PART(end,[1 3 5 11 13]), Y_PART(end,[1 3 5 11 13]+1), zeros(1,5), lb]';
            Y_PART(end,12) = Func_alphaB_VA_v2(Pv3); % velocity reset at touchdown
        end
        
        if iEVENT(i)==7
            % If this is EVENT 3, append right front touchdown FR:
            T_PART=[T_PART;T_PART(end)];
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv4 = [Y_PART(end,[1 3 5 11 13]), Y_PART(end,[1 3 5 11 13]+1), zeros(1,5), lb]';
            Y_PART(end,14) = Func_alphaF_VA_v2(Pv4); % velocity reset at touchdown
        end
        
        
        if iEVENT(i)==2
            % If this is EVENT 2, append left hind liftoff BL:
            T_PART=[T_PART;T_PART(end)]; 
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv1 = [Y_PART(end,[1 3 5 7 9]), Y_PART(end,[1 3 5 7 9]+1), zeros(1,5), lb]';
            Y_PART(end,8) = Func_alphaB_VA_v2(Pv1); % velocity reset at leftoff
        end
        
        if iEVENT(i)==4
            % If this is EVENT 4, append left front liftoff FL:
            T_PART=[T_PART;T_PART(end)];
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv2 = [Y_PART(end,[1 3 5 7 9]), Y_PART(end,[1 3 5 7 9]+1), zeros(1,5), lb]';
            Y_PART(end,10) = Func_alphaF_VA_v2(Pv2); % velocity reset at leftoff
        end
        
        if iEVENT(i)==6
            % If this is EVENT 6, append right hind liftoff BR:
            T_PART=[T_PART;T_PART(end)]; 
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv3 = [Y_PART(end,[1 3 5 11 13]), Y_PART(end,[1 3 5 11 13]+1), zeros(1,5), lb]';
            Y_PART(end,12) = Func_alphaB_VA_v2(Pv3); % velocity reset at leftoff
        end
        
        if iEVENT(i)==8
            % If this is EVENT 8, append right front liftoff FR:
            T_PART=[T_PART;T_PART(end)];
            Y_PART=[Y_PART;Y_PART(end,:)];
            Pv4 = [Y_PART(end,[1 3 5 11 13]), Y_PART(end,[1 3 5 11 13]+1), zeros(1,5), lb]';
            Y_PART(end,14) = Func_alphaF_VA_v2(Pv4); % velocity reset at leftoff
        end
        
    
        % Compose total solution
        T = [T;T_PART];
        Y = [Y;Y_PART];

        % Extract values at Events
        Y_EVENT(iEVENT(i),:)=Y(end,:);
        % Prepare initial values for next integration:
        T_START = T(end);
        Y_START = Y(end,:);
    end

    % Return the parameter set
    P = [tBL_TD_,tBL_LO_,tFL_TD_,tFL_LO_,tBR_TD_,tBR_LO_,tFR_TD_,tFR_LO,tAPEX,k,ks,J,l,osa,lb,kr];
    
    % Compute ground reation force
    [GRF,GRF_X,GRF_Y] = ComputeGRF(P,Y,T);
    GRFs = [GRF GRF_X GRF_Y];
     
    %**********************************************************************
    % Compute Residuals
    %**********************************************************************
    % Relbbel event values:
    YBL_TD = Y_EVENT(1,:)';
    YBL_LO = Y_EVENT(2,:)';
    YFL_TD = Y_EVENT(3,:)';
    YFL_LO = Y_EVENT(4,:)';

    YBR_TD = Y_EVENT(5,:)';
    YBR_LO = Y_EVENT(6,:)';
    YFR_TD = Y_EVENT(7,:)';
    YFR_LO = Y_EVENT(8,:)';    
    
    YAPEX  = Y_EVENT(9,:)';
    % Compute residuals
    residual = zeros(N+12,1);

    % Periodicity:
    residual(1:N-1) = Y(1,2:N).' - YAPEX(2:N); % x is not periodic
    
         
    % At the touch-down events, the feet have to be on the ground:
    residual(N+0)  = YBL_TD(3)  - lb*sin(YBL_TD(5)) - l*cos(YBL_TD(5)+YBL_TD(7));
    residual(N+1)  = YFL_TD(3)  + (1-lb)*sin(YFL_TD(5)) - l*cos(YFL_TD(5)+YFL_TD(9));
    % At the lift-off events, the feet also have to be on the ground:
    residual(N+2)  = YBL_LO(3)  - lb*sin(YBL_LO(5)) - l*cos(YBL_LO(5)+YBL_LO(7));
    residual(N+3)  = YFL_LO(3)  + (1-lb)*sin(YFL_LO(5)) - l*cos(YFL_LO(5)+YFL_LO(9));
    
    % At the touch-down events, the feet have to be on the ground:
    residual(N+4)  = YBR_TD(3)  - lb*sin(YBR_TD(5)) - l*cos(YBR_TD(5)+YBR_TD(11));
    residual(N+5)  = YFR_TD(3)  + (1-lb)*sin(YFR_TD(5)) - l*cos(YFR_TD(5)+YFR_TD(13));
    % At the lift-off events, the feet also have to be on the ground:
    residual(N+6)  = YBR_LO(3)  - lb*sin(YBR_LO(5)) - l*cos(YBR_LO(5)+YBR_LO(11));
    residual(N+7)  = YFR_LO(3)  + (1-lb)*sin(YFR_LO(5)) - l*cos(YFR_LO(5)+YFR_LO(13));
    
    %    Poincare section
    residual(end+1) = YAPEX(4);

%     % For infinite inertia
%     residual(end+1) = YAPEX(5)-0.00;
    
%     % Symmetry constraints(Front-Back)
%     residual(N+8)  = YBL_TD(7)  + YFL_LO(9);
%     residual(N+9)  = YFL_TD(9)  + YBL_LO(7);
%     residual(N+10) = YBL_TD(11) + YFL_LO(13);
%     residual(N+11) = YFL_TD(13) + YBL_LO(11);

%     % Symmetry constraints(Left-Right)
%     residual(end+1) = YFL_TD(9) + YFR_LO(13);
%     residual(end+1) = YFL_LO(9) + YFR_TD(13);
%     residual(end+1) = YBL_TD(7) + YBR_LO(11);
%     residual(end+1) = YBL_LO(7) + YBR_TD(11);

%     % Left and right symmetry:
%     residual(end+1)  = abs(YFL_TD(9)) - abs(YFR_TD(13));
%     residual(end+1)  = abs(YBL_TD(7)) - abs(YBR_TD(11));

%     % Front and back symmetry:
%     residual(end+1)  = YFL_TD(9) - YBL_TD(7);
%     residual(end+1)  = YFR_TD(13) - YBR_TD(11);
     
%     % Diagnal legs symmetry:
%     residual(end+1)  = YFL_TD(9) + YBR_LO(11);
%     residual(end+1)  = abs(YBL_TD(7)) - abs(YFR_TD(13));
    


    %**********************************************************************
    % Dynamics Function
    %**********************************************************************
    function dydt_ = ode(~,Y)
    % Extract individual states:
        x       = Y(1);
        dx      = Y(2);
        y       = Y(3);
        dy      = Y(4);
        phi     = Y(5);
        dphi    = Y(6);
        alphaBL  = Y(7);
        dalphaBL = Y(8);
        alphaFL  = Y(9);
        dalphaFL = Y(10);
        alphaBR  = Y(11);
        dalphaBR = Y(12);
        alphaFR  = Y(13);
        dalphaFR = Y(14);
        
        pos0 = [x;y];
        posB = pos0 + lb*[cos(phi + pi);sin(phi + pi)] ;
        posF = pos0 + (1-lb)*[cos(phi);sin(phi)] ;
         
        % Compute forces acting on the main body (only legs in contact
        % contribute): 

        BLforce = 0;
        FLforce = 0;
        BRforce = 0;
        FRforce = 0;
        % force applied on each leg, compute only in contact
        if contactBL
            BLforce = (l- posB(2)/cos(alphaBL+phi))*kb;
        end
        if contactFL
            FLforce = (l- posF(2)/cos(alphaFL+phi))*kf;
        end
        if contactBR
            BRforce = (l- posB(2)/cos(alphaBR+phi))*kb;
        end
        if contactFR
            FRforce = (l- posF(2)/cos(alphaFR+phi))*kf;
        end
        
        % Toral force and torque applied on the torso
        Fx  = -BLforce*sin(alphaBL+phi) - FLforce*sin(alphaFL+phi)...
              -BRforce*sin(alphaBR+phi) - FRforce*sin(alphaFR+phi);
        Fy  =  BLforce*cos(alphaBL+phi) + FLforce*cos(alphaFL+phi)...
              +BRforce*cos(alphaBR+phi) + FRforce*cos(alphaFR+phi);
        Tor = -BLforce*lb*cos(alphaBL)  + FLforce*(1-lb)*cos(alphaFL) ...
              -BRforce*lb*cos(alphaBR)  + FRforce*(1-lb)*cos(alphaFR);

        % Compute main body acceleration:
        ddx   = Fx;
        ddy   = Fy-1;
        ddphi = Tor/J;
       

        % Compute leg dynamics:

        AsvL = [x y phi alphaBL alphaFL dx dy dphi dalphaBL dalphaFL ddx ddy ddphi 0 0, lb]';  
    
        if contactBL % stance dynamics
            [dalphaBL,ddalphaBL] = Func_alphaB_VA_v2(AsvL);  
        else         % swing dynamics
            ddalphaBL = - ( Tor/J - Tor*lb*sin(alphaBL)/(J*l) ...
                         + Fx*cos(alphaBL + phi)/(M*l)...
                         + Fy*sin(alphaBL + phi)/(M*l)...
                         + alphaBL*ks/(l^2)... 
                         + dphi^2*lb*cos(alphaBL)/l);                       
        end
        
        if contactFL % stance dynamics
            [dalphaFL,ddalphaFL] = Func_alphaF_VA_v2(AsvL);
        else         % swing dynamics
            ddalphaFL = - ( Tor/J + Tor*(1-lb)*sin(alphaFL)/(J*l)...
                         + Fx*cos(alphaFL + phi)/(M*l)...
                         + Fy*sin(alphaFL + phi)/(M*l)...
                         + alphaFL*ks/(l^2)...
                         - dphi^2*(1-lb)*cos(alphaFL)/l);
        end


        AsvR = [x y phi alphaBR alphaFR dx dy dphi dalphaBR dalphaFR ddx ddy ddphi 0 0, lb]';  
     
        if contactBR % stance dynamics
            [dalphaBR,ddalphaBR] = Func_alphaB_VA_v2(AsvR);  
        else         % swing dynamics
            ddalphaBR = - ( Tor/J - Tor*lb*sin(alphaBR)/(J*l) ...
                         + Fx*cos(alphaBR + phi)/(M*l)...
                         + Fy*sin(alphaBR + phi)/(M*l)...
                         + alphaBR*ks/(l^2)... 
                         + dphi^2*lb*cos(alphaBR)/l);                       
        end
        
        if contactFR % stance dynamics
            [dalphaFR,ddalphaFR] = Func_alphaF_VA_v2(AsvR);
        else         % swing dynamics
            ddalphaFR = - ( Tor/J + Tor*(1-lb)*sin(alphaFR)/(J*l)...
                         + Fx*cos(alphaFR + phi)/(M*l)...
                         + Fy*sin(alphaFR + phi)/(M*l)...
                         + alphaFR*ks/(l^2)...
                         - dphi^2*(1-lb)*cos(alphaFR)/l);
        end
        
        
        dydt_ = [dx;ddx;dy;ddy;dphi;ddphi;...
                 dalphaBL;ddalphaBL;dalphaFL;ddalphaFL;...
                 dalphaBR;ddalphaBR;dalphaFR;ddalphaFR];
        
    end

end
%% Compute Stance Leg Dynamics
    function [dalphaB,ddalphaB] = Func_alphaB_VA_v2(in1)
    %Func_alphaB_VA_v2_V2
    %    [DALPHAB,DDALPHAB] = Func_alphaB_VA_v2_V2(IN1)

    %    This function was generated by the Symbolic Math Toolbox version 8.7.
    %    16-Feb-2023 22:36:05

    alphaB = in1(4,:);
    dalphaB = in1(9,:);
    ddphi = in1(13,:);
    ddx = in1(11,:);
    ddy = in1(12,:);
    dphi = in1(8,:);
    dx = in1(6,:);
    dy = in1(7,:);
    lb = in1(16,:);
    phi = in1(3,:);
    y = in1(2,:);
    t2 = cos(alphaB);
    t3 = sin(phi);
    t4 = alphaB+phi;
    t5 = alphaB.*2.0;
    t6 = alphaB.*3.0;
    t7 = dalphaB.^2;
    t8 = dphi.^2;
    t9 = phi.*2.0;
    t10 = phi.*3.0;
    t11 = cos(t4);
    t12 = sin(t4);
    t13 = phi+t4;
    t16 = t4.*2.0;
    dalphaB = -(dx+dphi.*y.*2.0+dx.*cos(t16)+dy.*sin(t16)-dphi.*lb.*t3-dphi.*lb.*sin(alphaB+t4))./(y.*2.0-lb.*t3.*2.0);
    if nargout > 1
        t18 = t4.*3.0;
        t14 = sin(t13);
        t15 = cos(t13);
        t17 = alphaB+t16;
        ddalphaB = -(ddx.*t11.*3.0+ddy.*t12+ddx.*cos(t18)+ddy.*sin(t18)+lb.*t8.*cos(t17)-ddphi.*lb.*sin(t17)+dalphaB.*dy.*t11.*8.0+dphi.*dy.*t11.*8.0-ddphi.*lb.*t14+ddphi.*t11.*y.*4.0-lb.*t2.*t7.*4.0-lb.*t2.*t8.*6.0+lb.*t7.*t15.*4.0+lb.*t8.*t15+t7.*t12.*y.*8.0+t8.*t12.*y.*8.0-dalphaB.*dphi.*lb.*t2.*1.2e+1+dalphaB.*dphi.*lb.*t15.*4.0+dalphaB.*dphi.*t12.*y.*1.6e+1)./(lb.*t14.*-2.0+t11.*y.*4.0+lb.*sin(alphaB).*2.0);
    end

    end
    function [dalphaF,ddalphaF] = Func_alphaF_VA_v2(in1)
    %FUNC_ALPHAF_VA_V2
    %    [DALPHAF,DDALPHAF] = FUNC_ALPHAF_VA_V2(IN1)

    %    This function was generated by the Symbolic Math Toolbox version 8.7.
    %    16-Feb-2023 22:36:08

    alphaF = in1(5,:);
    dalphaF = in1(10,:);
    ddphi = in1(13,:);
    ddx = in1(11,:);
    ddy = in1(12,:);
    dphi = in1(8,:);
    dx = in1(6,:);
    dy = in1(7,:);
    lb = in1(16,:);
    phi = in1(3,:);
    y = in1(2,:);
    t2 = cos(phi);
    t3 = sin(phi);
    t4 = alphaF+phi;
    t6 = lb-1.0;
    t5 = tan(t4);
    t8 = t3.*t6;
    t7 = t5.^2;
    t10 = -t8;
    t12 = t2.*t5.*t6;
    t16 = -1.0./(t8-y);
    t9 = t7+1.0;
    t14 = t10+y;
    t15 = -t12;
    t11 = dy.*t9;
    t13 = 1.0./t9;
    dalphaF = (t13.*(dx+dy.*t5-dphi.*(t10+t12+t9.*(t8-y))))./(t8-y);
    if nargout > 1
        t17 = -t9.*(t8-y);
        t18 = t5.*t9.*(t8-y).*-2.0;
        t19 = dalphaF.*t18;
        ddalphaF = (t13.*(ddx+dphi.*(t11+t19+dphi.*(t18+t2.*t6+t5.*t8-t2.*t6.*t9.*2.0)-dalphaF.*t2.*t6.*t9)-dalphaF.*(-t11+dphi.*(t2.*t6.*t9+t5.*t9.*(t8-y).*2.0)+dalphaF.*t5.*t9.*(t8-y).*2.0)+ddy.*t5-ddphi.*(t10+t12+t9.*(t8-y))+dy.*(dalphaF.*t9+dphi.*t9)))./(t8-y);
    end

    end
%%  Compute Ground Reaction Force
function [GRF,GRF_X,GRF_Y] = ComputeGRF(P,Y,T)

    tBL_TD = P(1);
    tBL_LO = P(2);
    tFL_TD = P(3);
    tFL_LO = P(4);
    tBR_TD = P(5);
    tBR_LO = P(6);
    tFR_TD = P(7);
    tFR_LO = P(8);
    k  = P(10);
    l  = P(13);
    lb = P(15);
    kr = P(16);
    
    kb = 2*k/(1+ 1/kr);
    kf = 2*k/(1+ kr);
    
    n = length(T);
    
    FBLx = zeros(n,1);
    FFLx = zeros(n,1);
    FBRx = zeros(n,1);
    FFRx = zeros(n,1);
    
    FBLy = zeros(n,1);
    FFLy = zeros(n,1);
    FBRy = zeros(n,1);
    FFRy = zeros(n,1);
    
    FBL = zeros(n,1);
    FFL = zeros(n,1);
    FBR = zeros(n,1);
    FFR = zeros(n,1);
    
    
    % Compute vertical GRF using the data   
    for  i = 1:n

        t_ = T(i);
        if ((t_>tBL_TD && t_<tBL_LO && tBL_TD<tBL_LO) || ((t_<tBL_LO || t_>tBL_TD) && tBL_TD>tBL_LO))
            contactBL = true;
        else
            contactBL = false;
        end
        if ((t_>tFL_TD && t_<tFL_LO && tFL_TD<tFL_LO) || ((t_<tFL_LO || t_>tFL_TD) && tFL_TD>tFL_LO))
            contactFL = true;
        else
            contactFL = false;
        end


        if ((t_>tBR_TD && t_<tBR_LO && tBR_TD<tBR_LO) || ((t_<tBR_LO || t_>tBR_TD) && tBR_TD>tBR_LO))
            contactBR = true;
        else
            contactBR = false;
        end
        if ((t_>tFR_TD && t_<tFR_LO && tFR_TD<tFR_LO) || ((t_<tFR_LO || t_>tFR_TD) && tFR_TD>tFR_LO))
            contactFR = true;
        else
            contactFR = false;
        end  

        y_ = Y(i,:);
        
        x        = y_(1);
        y        = y_(3);
        phi      = y_(5);        
        alphaBL   = y_(7);
        alphaFL   = y_(9);
        alphaBR   = y_(11);
        alphaFR   = y_(13);


        pos0 = [x;y];
        posB = pos0 + lb*[cos(phi + pi);sin(phi + pi)] ;
        posF = pos0 + (1-lb)*[cos(phi);sin(phi)] ;

        % Compute forces acting on the main body (only legs in contact
        % contribute): 
        
        if contactBL
            FBLx(i) = (l- posB(2)/cos(alphaBL+phi))*kb*sin(alphaBL+phi);
        end
        if contactFL
            FFLx(i) = (l- posF(2)/cos(alphaFL+phi))*kf*sin(alphaFL+phi);
        end 
        if contactBR
            FBRx(i) = (l- posB(2)/cos(alphaBR+phi))*kb*sin(alphaBR+phi);
        end
        if contactFR
            FFRx(i) = (l- posF(2)/cos(alphaFR+phi))*kf*sin(alphaFR+phi);
        end
              
        if contactBL
            FBLy(i) = (l- posB(2)/cos(alphaBL+phi))*kb*cos(alphaBL+phi);
        end
        if contactFL
            FFLy(i) = (l- posF(2)/cos(alphaFL+phi))*kf*cos(alphaFL+phi);
        end 
        if contactBR
            FBRy(i) = (l- posB(2)/cos(alphaBR+phi))*kb*cos(alphaBR+phi);
        end
        if contactFR
            FFRy(i) = (l- posF(2)/cos(alphaFR+phi))*kf*cos(alphaFR+phi);
        end
        
        if contactBL
            FBL(i) = (l- posB(2)/cos(alphaBL+phi))*kb;
        end
        if contactFL
            FFL(i) = (l- posF(2)/cos(alphaFL+phi))*kf;
        end 
        if contactBR
            FBR(i) = (l- posB(2)/cos(alphaBR+phi))*kb;
        end
        if contactFR
            FFR(i) = (l- posF(2)/cos(alphaFR+phi))*kf;
        end   

    end
    GRF_X = [FBLx FFLx FBRx FFRx];
    GRF_Y = [FBLy FFLy FBRy FFRy];
    GRF   = [FBL  FFL  FBR  FFR ];
end
