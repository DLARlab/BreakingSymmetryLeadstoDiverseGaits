%% Identifying the type of gaits for a given solution or solution branch, return the name of gait and the color for visualization
function [gait,abbr, color_plot, linetype] = Gait_Identification(results)
    if size(results,2)==1
        [gait,abbr, color_plot, linetype] = Type_of_Gait(results);
    else
        downsample_rate = fix(size(results,2)/10);
        indices = downsample(1:size(results,2),downsample_rate);
        count = 1;
        for i = 1:length(indices)
            X = results(:,indices(i));
            [gait_current,abbr_current, color_current, linetype_current] = Type_of_Gait(X);
            if i>2
                if string(gait_current)==string(gait_last)
                    count = count +1;
                    if count>5
                        gait = gait_current;
                        abbr = abbr_current;
                        color_plot = color_current;
                        linetype = linetype_current;
                        break
                    elseif i > (length(indices)-1)
                        disp('Gait identification error for current branch.')
                    end
                end
            end
            gait_last = gait_current;
        end

    end
end

%% Gait identifying function
function [gait,abbr, color_plot, linetype] = Type_of_Gait(X)

    threshold = 1e-6; % threshold for identifying 'equal'

    if abs(X(14)-X(18))<threshold && abs(X(16)-X(20))<threshold && abs(X(14)-X(16))<threshold && abs(X(16)-X(20))<threshold
        type1 = 'Pronking';
        abbr1 = 'PF';
        color_plot = [0 0.4470 0.7410];
    elseif abs(X(14)-X(18))<threshold && abs(X(16)-X(20))<threshold
        type1 = 'Bounding';
        abbr1 = 'B';
        color_plot = [0.8500 0.3270 0.0980];
    elseif abs(X(14)-X(18))<threshold 
        type1 = 'Half-Bounding with Front Legs Spread';
        abbr1 = 'F';
        color_plot = [0.4660 0.6740 0.1880];
    elseif abs(X(16)-X(20))<threshold 
        type1 = 'Half-Bounding with Hind Legs Spread';
        abbr1 = 'H';
        color_plot = [0.9290 0.6940 0.1270];
    else
        type1 = 'Galloping';
        abbr1 = 'G';
        color_plot = [0.4940 0.1840 0.5560];
    end
    
    if abs(X(5))<1e-9
        type2 = '';
        abbr2 = '';
        linetype = '-';
    elseif X(5)>0
        type2 = '_with Gathered Suspension';
        abbr2 = 'G';
        linetype = '-';
    elseif  X(5)<0
        type2 = '_with Extended Suspension';
        abbr2 = 'E';
        linetype = '--';
    end

    gait = string(type1) + string(type2);
    abbr = string(abbr1) + string(abbr2);

end