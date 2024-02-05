function [ x1, y1, f,v ] = ComputeBodyGraphics(BodyJPos,lb)

l  = 1.2; % length of the rectangle
h  = 0.4; % width of the rectangle
n  = 30;  % number of point along length
m  = 10;  % number of point along width

if lb == 0.5
    % % Create vertices on the rectangle
    v_ = [ repmat([-l/2, h/2],n,1) + [linspace(0,l,n)',zeros(n,1)];                                                                             
           repmat([ l/2, h/2],m,1) - [zeros(m,1),linspace(0,h,m)'];                                                                                           
           repmat([ l/2,-h/2],n,1) - [linspace(0,l,n)',zeros(n,1)];                                                                                                                      
           repmat([-l/2,-h/2],m,1) + [zeros(m,1),linspace(0,h,m)']; ];
else
    h  = 0.45;
    offset = sign(0.5-lb)*sqrt(abs((0.5-lb)*10))/15;
    v_ = [ repmat([-l/2, h/2+offset],n,1) + [linspace(0,l,n)',linspace(0,-2*offset,n)'];                                                                             
           repmat([ l/2, h/2-offset],m,1) - [zeros(m,1),linspace(0,h-2*offset,m)'];                                                                                           
           repmat([ l/2,-h/2+offset],n,1) - [linspace(0,l,n)',linspace(0,+2*offset,n)'];                                                                                                                      
           repmat([-l/2,-h/2-offset],m,1) + [zeros(m,1),linspace(0,h+2*offset,m)']; ];

    v_ (:,1)= v_ (:,1) + (0.5-lb);
end

% Shift center to x,y and rotate for phi:
v  = LTrans(v_',BodyJPos(3),BodyJPos(1:2)');
% Manually find points to create stripes
f = [ [4:3:n,                n+2:2:n+10-2    ]; 
      [2*m+2*n-2:-2:m+2*n+2, m+2*n-2:-3:m+n+2] ]';
% Four vertices of the rectangle
x1 = [v(1,1),v(1+n,1),v(1+n+m,1),v(1+2*n+m,1),v(1,1)]; 
y1 = [v(1,2),v(1+n,2),v(1+n+m,2),v(1+2*n+m,2),v(1,2)];

% *********************************************************************
% Local function, rotate and shift objects using linear transformation
function VecTrans = LTrans(VecTrans, gamma, vecS)
    % Get absolute angle of each object and the vector of translational
    % motion. VecTrans is a 2 by n vector and the output VecTrans is a
    % n by 2 vector
    % (which must be on the MATLAB search path):
    % Rotational matrix:
    RotM = [ cos(gamma), -sin(gamma);
             sin(gamma),  cos(gamma)];
    vec_rot = RotM*VecTrans;    % Rotate 
    vec_ = vec_rot + vecS*ones(1,size(VecTrans,2)); % Shift
    VecTrans = vec_'; % Change to column vector.
end
% *********************************************************************
end

