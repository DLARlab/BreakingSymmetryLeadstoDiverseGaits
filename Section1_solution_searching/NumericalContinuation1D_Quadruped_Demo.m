function [results,flag] = NumericalContinuation1D_Quadruped_Demo(X1,X2,Para,radius,numOPTS)
    % Define the solution matrix. Search negative direction first.
    results(:,1) = [X1;Para]; 
    results(:,2) = [X2;Para];
    % results = flip(sortrows(results.',1).',2); % sort so that result2 < result1: searching for low speed first
    results = sortrows(results.',1).';      % sort so that result1 < result2: searching for high speed first

    % Parameter define temporary storage or not
    SaveTempSol = 0;
    % Illustrate the process of searching 
    IlluSols = 1;
    
    if IlluSols==1
        BranchPlot = BranchPlotSearching(results);
    end
     
     % Searching for solutions along 1 direction
     [results,flag] = UniDirectionSearch(results,radius,numOPTS,SaveTempSol,IlluSols,BranchPlot);

                
end
%% Unidirection Search
function [results,flag] = UniDirectionSearch(results, Radius, numOPTS, SaveTempSol, IlluSols, BranchPlot)
    x_Last = results(1:22,end);
    x_Current = x_Last;
    Para = results(23:end,1);
    ns = size(results,2);
for k = 1:20
       % Set special options for the stiff range
%        if x_Current(1)>15
%            OPTS = optimset('Algorithm','levenberg-marquardt',... 
%                    'Display','iter',...
%                    'MaxFunEvals',50000,...
%                    'MaxIter',3000,...
%                    'UseParallel', false,...
%                    'TolFun',1e-12,...
%                    'TolX',1e-12);
%        else % Use the original settings
           OPTS = numOPTS;         
%        end
       % Numerical Continuation using Prediction-Correction method: assume a new solution always exists in the neighborhood of existing solutions
       % Prediction of new solution: X(k)* = X(k-1) + (X(k-1)-X(k-2))*d
       % Prediction of the k th solution
       X0 = x_Current + (results(1:22,k+ns-1)-results(1:22,k+ns-2));
       
       if IlluSols == 1
           figure(BranchPlot.fig)
           current_speed_range = [BranchPlot.axes.XLim(1) X0(1) results(1,end)];
           if sign(min(current_speed_range)) == 1 
                BranchPlot.axes.XLim = [min(current_speed_range)  max(current_speed_range)*(1+sign(max(current_speed_range))*0.1)];
           else
                BranchPlot.axes.XLim = [min(current_speed_range)*(1-sign(min(current_speed_range))*0.1)  BranchPlot.axes.XLim(2)];
           end
           if k == 1
               predict_solution = scatter3(X0(1),X0(3),X0(5), 120, 'Parent',BranchPlot.axes,'MarkerEdgeColor',[0 0 0]);
               predict_text_position =  BranchPlot.Text.textPO.Position + [results(1,k+ns-1)-results(1,k+ns-2) 0  results(5,k+ns-1)-results(5,k+ns-2)+0.125*BranchPlot.axes.ZLim(2)];
               predict_text =  text( predict_text_position(1),predict_text_position(2),predict_text_position(3),...
                                     'Predicted Solution','HorizontalAlignment','center','Parent',BranchPlot.axes);
           else
               [gait,abbr, color_plot, linetype] = Gait_Identification(xFINAL);
               scatter3(results(1,end),results(3,end),results(5,end),120,'filled','Parent',BranchPlot.axes,...
                                           'MarkerEdgeColor',color_plot,'MarkerFaceColor',color_plot)
               BranchPlot.Text.textPO.Position(1) = results(1,end);
               predict_solution.XData = predict_solution.XData + results(1,k+ns-1)-results(1,k+ns-2);
               predict_solution.ZData = predict_solution.ZData + results(5,k+ns-1)-results(5,k+ns-2);
               predict_text.Position = predict_text.Position + [results(1,k+ns-1)-results(1,k+ns-2) 0  results(5,k+ns-1)-results(5,k+ns-2)];
           end
       end
       % Correction of the prediction of k th solution
       [xFINAL,fval] = ContinuationFun_Quadruped_v2(X0,x_Last,Para,Radius,OPTS);
       
              
       % In any other cases, set the flag of ternimation to be:
       flag = 'Algorithm reach the maximum number of iterations.';
       
       % Save the solution and keep searching if the solution is periodic,
       % or stop searching when the solution is not converging
       if ( norm(fval) < 1e-9 && 0<xFINAL(1)<15 )|| ( norm(fval) < 1e-6 && xFINAL(1)>15 ) % && xFINAL(1)<18) || (norm(fval) < 5e-7 && xFINAL(1)>18)
          % The solution is periodic
          x_Current = xFINAL;  
          x_Last = x_Current;
          results(:,k+ns) = [x_Current;Para];
          
          % Save the solution and continue searching
          if SaveTempSol == 1
              save('solution_staging.mat','results')
          end
          % Display the current solution
          pause(0.5); fprintf('\n') ; 
          pause(0.5); disp('Current Speed = '); disp(x_Current(1));
          pause(0.5); disp('Residual = ');      disp(norm(fval))
          pause(0.5)
       elseif ~(string(TimmingBoundaryCheck(results(1:22,end)))=='Null')
          % Terminated because the current solution is reaching the timming boundary.
          flag = 'Algorithm stops because event timing reaching the ' + string(TimmingBoundaryCheck(results(1:22,end))) +' of the simulation.';         
          pause(3)
          break
       else % Terminated because the current solution cannot converge due to the numerical issue based the ode options.
          flag = 'Algorithm stops because of numerical issues.';
          pause(3)
          break;
       end
       
                  
       % There are additional criterias that will terminate the algorithm:
       % 1.Solution that have horizontal speed < 0 has no physical meaning    
       if xFINAL(1)< -0.01
            flag = 'Algorithm stopped: Torso speed smaller than zero.';
            break;
       end
       % 2. Solution branch shows in a circle
       dir = sign(results(1,2)-results(1,1));
       if sign(xFINAL(1) - results(1,k+ns-1))==dir && norm(xFINAL(1:13)-results(1:13,1)) < Radius  && dir*xFINAL(1) > dir*results(1,1)
            flag = 'Algorithm stopped: Solution Circle Found.';
            break;
       end
       % 3. Reaching a bifurcation point: pairs of Couplet increases
       if Couplet(xFINAL)>Couplet(results(1:22,1))
            flag = 'Algorithm stopped: Reaching a Bifurcation Point.';
            break;
       end

       % 4.Solution that have horizontal speed > 15 has no physical meaning    
       if xFINAL(1)> 10
            flag = 'Algorithm stopped: Torso speed faster than 15.';
            disp('Algorithm stopped: Torso speed faster than 15.')
            break;
       end
       
end
disp(flag)
end

%% Continuation Function that is used to set the 1-D continuation
function [xFINAL,fval] = ContinuationFun_Quadruped_v2(X0,x_Last,Para,radius,numOPTS)

% Find solution X zero in residual values
[xFINAL,fval] = fsolve(@Residual, X0, numOPTS);


function residual = Residual(X)
    
    % Build new residual values zero in fsolve
    [residual,T,Y,P] = Quadrupedal_ZeroFun_v2(X,Para);
    % Constrain that limit the distance of new solution from current one
    residual(end+1) = norm(X(1:13)-x_Last(1:13))-radius;
end

end

%% Couplet check of leg pairs: Synchronization between left and right legs
function noc = Couplet(X)
noc = 0;
if norm(X(14)-X(18))<0.001 && norm(X(15)-X(19))<0.001
     noc = noc + 1;
end
if norm(X(16)-X(20))<0.001 && norm(X(17)-X(21))<0.001
     noc = noc + 1;
end
end

%% Plot animation for the solution searching process
function BranchPlot = BranchPlotSearching(results)
       % Plot position settings
        ScreenSize = get(0,'ScreenSize');
        PlotSize = [(2.5/10)*ScreenSize(3)   (12/16)*(2.5/10)*ScreenSize(3)
                 ScreenSize(3)/5        (12/16)*ScreenSize(3)/5     ];
        PlotPositions = [(1/2)*ScreenSize(3)-(2/2 + 1/10)*PlotSize(1,1)  (0.5/10)*ScreenSize(4)+PlotSize(2,2)  PlotSize(1,1)  PlotSize(1,2)
                      (1/2)*ScreenSize(3)-(0/2 - 1/10)*PlotSize(1,1)  (0.5/10)*ScreenSize(4)+PlotSize(2,2)  PlotSize(1,1)  PlotSize(1,2)];
        % Simulate the sytem
        
        [residual,T,Y,P,GRF,Y_EVENT] = Quadrupedal_ZeroFun_v2(results(1:22,end),results(23:end,end));
        [gait,abbr, color_plot, linetype] = Gait_Identification(results(1:22,end));
        figure(206)
        BranchPlot =  SLIP_PeriodicOrbit_Quad(Y,PlotPositions(2,:),figure(206),color_plot);
        
        % Update the current state on the orbit
        n = round(T(end)*100); % # of frames per step
        tFrame = linspace(0, T(end), n+1);
        for j = 1:n
            % set y date for updating animation and periodic orbit 
            y    = interp1(T' + linspace(0,1e-5,length(T)), Y,   tFrame(j));   
            BranchPlot.update(y);
            pause(0.02)
        end
        % Change the view
        figure(BranchPlot.fig)
        [caz,cel] = view;
        for j = 1:n
            dview = 0:90/n:90;
            % set y date for updating animation and periodic orbit 
            y    = interp1(T' + linspace(0,1e-5,length(T)), Y,   tFrame(j));   
            BranchPlot.update(y);
            view([caz-dview(j)/2 cel-dview(j)/2])
            pause(0.02)
        end
        
        pause(2)
        delete(BranchPlot.Orbit)
        BranchPlot.Text.textPO.Position(1) = 0;
        BranchPlot.Text.textPO.Position = BranchPlot.Text.textPO.Position + [results(1,end) 0  0.25*BranchPlot.axes.ZLim(2)];
        BranchPlot.Text.textPO.String = 'Current Solution';
        pause(2)
        delete(BranchPlot.Poincare_Section)
        delete(BranchPlot.Text.textPS)
        
        current_speed_range = [BranchPlot.axes.XLim(1) BranchPlot.axes.XLim(2) results(1,end-1) results(1,end)];   
        BranchPlot.axes.XLim = [min(current_speed_range*(1-sign(min(current_speed_range))*0.1))   max(current_speed_range)*(1+sign(max(current_speed_range)*0.1))];
        BranchPlot.axes.Title.String = 'Searching for Solutions on the Branch';
        scatter3(results(1,end-1),results(3,end-1),results(5,end-1),120,'filled','Parent',BranchPlot.axes,...
                                           'MarkerEdgeColor',color_plot,'MarkerFaceColor',color_plot)
end
