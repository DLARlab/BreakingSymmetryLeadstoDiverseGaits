function [ LegLength, LegAngle, BodyJPos, BackJPos, FrontJPos] = ComputeJoint_LegLA(y,P,T)

% Get a mapping for the state and parameter vectors.  This allows us
% to use a more readable syntax: "y(contStateIndices.dy)" instead of
% "y(3)" while still operating with vectors and not with structs.
% We keep the index-structs in memory to speed up processing
% persistent contStateIndices  systParamIndices discStateIndices
% if isempty(contStateIndices) || isempty(discStateIndices) || isempty(systParamIndices) 
%     [~, ~, contStateIndices] = ContStateDefinition();
%     [~, ~, discStateIndices] = DiscStateDefinition();
%     [~, ~, systParamIndices] = SystParamDefinition();
% end

lb       = P(15); % half of the body length

l_l      = 1;
l_r      = 1;


tBL_TD = P(1);
tBL_LO = P(2);
tFL_TD = P(3);
tFL_LO = P(4);

tBR_TD = P(5);
tBR_LO = P(6);
tFR_TD = P(7);
tFR_LO = P(8);

tAPEX = P(9);

    if tBL_TD < 0
         tBL_TD = tBL_TD + tAPEX;
    end
    if tBL_TD > tAPEX
         tBL_TD = tBL_TD - tAPEX;
    end 
    
    if tBL_LO < 0
         tBL_LO = tBL_LO + tAPEX;
    end
    if tBL_LO > tAPEX
         tBL_LO = tBL_LO - tAPEX;
    end 
    
    if tFL_LO < 0
         tFL_LO = tFL_LO + tAPEX;
    end
    if tFL_LO > tAPEX
         tFL_LO = tFL_LO - tAPEX;
    end
    
    if tFL_TD < 0
         tFL_TD = tFL_TD + tAPEX;
    end
    if tFL_TD > tAPEX
         tFL_TD = tFL_TD - tAPEX;
    end
    
    
    
    if tBR_TD < 0
         tBR_TD = tBR_TD + tAPEX;
    end
    if tBR_TD > tAPEX
         tBR_TD = tBR_TD - tAPEX;
    end 
    
    if tBR_LO < 0
         tBR_LO = tBR_LO + tAPEX;
    end
    if tBR_LO > tAPEX
         tBR_LO = tBR_LO - tAPEX;
    end 
    
    if tFR_LO < 0
         tFR_LO = tFR_LO + tAPEX;
    end
    if tFR_LO > tAPEX
         tFR_LO = tFR_LO - tAPEX;
    end
    
    if tFR_TD < 0
         tFR_TD = tFR_TD + tAPEX;
    end
    if tFR_TD > tAPEX
         tFR_TD = tFR_TD - tAPEX;
    end
    
    
if ((T>tBL_TD && T<tBL_LO && tBL_TD<tBL_LO) || ((T<tBL_LO || T>tBL_TD) && tBL_TD>tBL_LO))
    contactBL = true;
else
    contactBL = false;
end
if ((T>tFL_TD && T<tFL_LO && tFL_TD<tFL_LO) || ((T<tFL_LO || T>tFL_TD) && tFL_TD>tFL_LO))
    contactFL = true;
else
    contactFL = false;
end

if ((T>tBR_TD && T<tBR_LO && tBR_TD<tBR_LO) || ((T<tBR_LO || T>tBR_TD) && tBR_TD>tBR_LO))
    contactBR = true;
else
    contactBR = false;
end
if ((T>tFR_TD && T<tFR_LO && tFR_TD<tFR_LO) || ((T<tFR_LO || T>tFR_TD) && tFR_TD>tFR_LO))
    contactFR = true;
else
    contactFR = false;
end


yx       = y(1);
yy       = y(3);
phi      = y(5);
alphaBL   = y(7);
alphaFL   = y(9);

alphaBR   = y(11);
alphaFR   = y(13);


x_b = -lb*cos(phi) + yx;
x_f = (1-lb)*cos(phi) + yx;
y_b = -lb*sin(phi) + yy;
y_f = (1-lb)*sin(phi) + yy;


if contactBL == true
    l_leg_BL = y_b/cos(alphaBL+phi);    
else    
    l_leg_BL = l_l;
end
gamma_leg_BL = alphaBL + phi; 

if contactFL == true
    l_leg_FL = y_f/cos(alphaFL+phi);    
else    
    l_leg_FL = l_r;
end
gamma_leg_FL = alphaFL + phi; 

if contactBR == true
    l_leg_BR = y_b/cos(alphaBR+phi);    
else    
    l_leg_BR = l_l;
end
gamma_leg_BR = alphaBR + phi; 

if contactFR == true
    l_leg_FR = y_f/cos(alphaFR+phi);    
else    
    l_leg_FR = l_r;
end
gamma_leg_FR = alphaFR + phi;                         


LegLength = struct('BL',l_leg_BL,'FL',l_leg_FL,'BR',l_leg_BR,'FR',l_leg_FR);
LegAngle  = struct('BL',gamma_leg_BL,'FL',gamma_leg_FL,'BR',gamma_leg_BR,'FR',gamma_leg_FR);
BodyJPos  = [yx, yy, phi];
BackJPos  = [x_b;y_b];
FrontJPos = [x_f;y_f];

end

