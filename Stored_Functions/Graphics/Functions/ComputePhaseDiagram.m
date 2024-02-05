function [vertices_box,pos_text,pos_PhaseBars] = ComputePhaseDiagram(x,P)
% A box will be drawn to include the phase bars and discriptions
% Define the center and size of box
x_box = x + 0.93;
y_box = 1.2+0.55;
l_box = 0.45;
w_box = 0.16;
l_text = 0.08;
% Define vertices of box 
box_x = [x_box-l_box-l_text,    x_box+l_box+0.5*l_text... 
         x_box+l_box+0.5*l_text,    x_box-l_box-l_text];
box_y = [y_box-w_box, y_box-w_box, y_box+w_box,  y_box+w_box];
vertices_box = [box_x; box_y];


% Define the origin, length, width, color of the phase bar
x_origin = x_box - l_box + 1.2*l_text;
len_whole = (l_box-0.05)*2;

width = 0.037;
yshift = 0.030;

% Show the descriptions for the phase bars
pos_text = [x_box-l_box-0.6*l_text,  y_box+1.5*yshift+1.5*width
            x_box-l_box-0.6*l_text,  y_box+0.5*yshift+0.5*width
            x_box-l_box-0.6*l_text,  y_box-0.5*yshift-0.5*width
            x_box-l_box-0.6*l_text,  y_box-1.5*yshift-1.5*width];
        



% Drawing the phase bars: draw two layers, one for stance face and one for flight phase.
% If tTD<tLO, draw flight phase(white) and then cover it with stance phase(black);
% If tTD>tLO, draw stance phase(white) and then cover it with flight phase(black).

% Define the color of phase bar
cp = [0 0 0];
% X = [0; 0; 0; 0; 0; 0; 0; tEvents];
% switch  Gaitidentify(X)
%     case 0
%         cp = [0.4940 0.1840 0.5560];
%     case 1
%         cp = [0.8500 0.3250 0.0980];
%     case 2
%         cp = [0 0.4470 0.7410];
%     case 3
%         cp = [0.4660 0.6740 0.1880];
%     case 4
%         cp = [0.9290 0.6940 0.1250];
% end

% Event timings that determine the sequence and the length of phase bars
tEvents = P(1:9);
% Back left phase bar
% vertices of base bar 
bar_blx_ =[x_origin,  x_origin+len_whole,   x_origin+len_whole,   x_origin];
bar_bly_ = [y_box + width + 1.5*yshift,     y_box + width + 1.5*yshift...
           y_box + 2*width + 1.5*yshift,    y_box + 2*width + 1.5*yshift];
if tEvents(1)<tEvents(2)
    % vertics of duration bar
    bar_blx = [x_origin + (tEvents(1)/tEvents(end))*len_whole,   x_origin + (tEvents(2)/tEvents(end))*len_whole...
              x_origin + (tEvents(2)/tEvents(end))*len_whole,    x_origin + (tEvents(1)/tEvents(end))*len_whole];
    bar_bly = [y_box + width + 1.5*yshift,      y_box + width + 1.5*yshift...
               y_box + 2*width + 1.5*yshift,    y_box + 2*width + 1.5*yshift];
elseif tEvents(1)>tEvents(2)
    bar_blx = [x_origin + (tEvents(2)/tEvents(end))*len_whole,   x_origin + (tEvents(1)/tEvents(end))*len_whole...
              x_origin + (tEvents(1)/tEvents(end))*len_whole,    x_origin + (tEvents(2)/tEvents(end))*len_whole];
    bar_bly = [y_box + width + 1.5*yshift,      y_box + width + 1.5*yshift...
               y_box + 2*width + 1.5*yshift,    y_box + 2*width + 1.5*yshift];
end   
pos_BL = [bar_blx_; bar_bly_; bar_blx; bar_bly];

% Back right phase bar
% vertices of base bar 
bar_brx_ = [x_origin,   x_origin+len_whole,   x_origin+len_whole,   x_origin];
bar_bry_ = [y_box - width - 1.5*yshift,       y_box - width - 1.5*yshift...
           y_box - 2*width - 1.5*yshift,      y_box - 2*width - 1.5*yshift];
if tEvents(5)<tEvents(6)
    % vertics of duration bar
    bar_brx = [x_origin + (tEvents(5)/tEvents(end))*len_whole,   x_origin + (tEvents(6)/tEvents(end))*len_whole...
              x_origin + (tEvents(6)/tEvents(end))*len_whole,    x_origin + (tEvents(5)/tEvents(end))*len_whole];
    bar_bry = [y_box - width - 1.5*yshift,         y_box - width - 1.5*yshift...
               y_box - 2*width - 1.5*yshift,       y_box - 2*width - 1.5*yshift];
elseif tEvents(5)>tEvents(6)
    bar_brx = [x_origin + (tEvents(6)/tEvents(end))*len_whole  x_origin + (tEvents(5)/tEvents(end))*len_whole...
              x_origin + (tEvents(5)/tEvents(end))*len_whole  x_origin + (tEvents(6)/tEvents(end))*len_whole];
    bar_bry = [y_box - width - 1.5*yshift,         y_box - width - 1.5*yshift...
               y_box - 2*width - 1.5*yshift,       y_box - 2*width - 1.5*yshift];
end
pos_BR = [bar_brx_; bar_bry_; bar_brx; bar_bry];


% Front left phase bar
% vertices of base bar 
bar_flx_ = [x_origin,   x_origin+len_whole,    x_origin+len_whole,    x_origin];
bar_fly_ = [y_box + 0*width + 0.5*yshift,      y_box + 0*width + 0.5*yshift...
            y_box + 1*width + 0.5*yshift,      y_box + 1*width + 0.5*yshift];
if tEvents(3)<tEvents(4)
    % vertics of duration bar    
    bar_flx = [x_origin + (tEvents(3)/tEvents(end))*len_whole,  x_origin + (tEvents(4)/tEvents(end))*len_whole...
              x_origin + (tEvents(4)/tEvents(end))*len_whole,   x_origin + (tEvents(3)/tEvents(end))*len_whole];
    bar_fly = [y_box + 0*width + 0.5*yshift,    y_box + 0*width + 0.5*yshift...
               y_box + 1*width + 0.5*yshift,    y_box + 1*width + 0.5*yshift];
elseif tEvents(3)>tEvents(4)   
    bar_flx = [x_origin + (tEvents(4)/tEvents(end))*len_whole,  x_origin + (tEvents(3)/tEvents(end))*len_whole...
              x_origin + (tEvents(3)/tEvents(end))*len_whole,   x_origin + (tEvents(4)/tEvents(end))*len_whole];
    bar_fly = [y_box + 0*width + 0.5*yshift,    y_box + 0*width + 0.5*yshift...
               y_box + 1*width + 0.5*yshift,    y_box + 1*width + 0.5*yshift];
end
pos_FL = [bar_flx_; bar_fly_; bar_flx; bar_fly];

% Front right phase bar
% vertices of base bar 
bar_frx_ = [x_origin,   x_origin+len_whole,    x_origin+len_whole,    x_origin];
bar_fry_ = [y_box - 0*width - 0.5*yshift,      y_box - 0*width - 0.5*yshift...
           y_box - 1*width - 0.5*yshift,       y_box - 1*width - 0.5*yshift];
if tEvents(7)<tEvents(8)
    % vertics of duration bar
    bar_frx = [x_origin + (tEvents(7)/tEvents(end))*len_whole,   x_origin + (tEvents(8)/tEvents(end))*len_whole...
              x_origin + (tEvents(8)/tEvents(end))*len_whole,    x_origin + (tEvents(7)/tEvents(end))*len_whole];
    bar_fry = [y_box - 0*width - 0.5*yshift,   y_box - 0*width - 0.5*yshift...
               y_box - 1*width - 0.5*yshift,   y_box - 1*width - 0.5*yshift];
elseif tEvents(7)>tEvents(8)
 
    bar_frx = [x_origin + (tEvents(8)/tEvents(end))*len_whole,   x_origin + (tEvents(7)/tEvents(end))*len_whole...
              x_origin + (tEvents(7)/tEvents(end))*len_whole,    x_origin + (tEvents(8)/tEvents(end))*len_whole];
    bar_fry = [y_box - 0*width - 0.5*yshift,    y_box - 0*width - 0.5*yshift...
               y_box - 1*width - 0.5*yshift,    y_box - 1*width - 0.5*yshift];
end
pos_FR = [bar_frx_; bar_fry_; bar_frx; bar_fry];

pos_PhaseBars = [pos_BL; pos_BR; pos_FL; pos_FR];

end

