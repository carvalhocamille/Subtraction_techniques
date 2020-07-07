%% BoundaryCurve.m
% Computes the parametric features of the boundary curve
% License MIT (see Readme.md)

function [ r, tangent, normal, Jacobian, kappa ] = BoundaryCurve( t )

% compute the parametric curve

% x1   = cos(t) + 0.65 * cos(2*t) - 0.65;
% x1p1 = -sin(t) - 1.3 * sin(2*t);
% x1p2 = -cos(t) - 2.6 * cos(2*t);
% 
% x2   = 1.5 * sin(t);
% x2p1 = 1.5 * cos(t);
% x2p2 = -1.5 * sin(t);
% 
% r   = [ x1   x2   ];
% rp1 = [ x1p1 x2p1 ];
% rp2 = [ x1p2 x2p2 ];
% 
% a = 3;
% b = 2;
% 
% r   = [  a * cos(t)  b * sin(t) ];
% rp1 = [ -a * sin(t)  b * cos(t) ];
% rp2 = [ -a * cos(t) -b * sin(t) ];

a = 0.4;
w = 5;

f   = 1.55 + a * cos( w * t );
fp1 = - w * a * sin( w * t );
fp2 = - w * w * a * cos( w * t );

r    = [ f .* cos(t) ...
    f .* sin(t) ];

rp1  = [ fp1 .* cos(t) - f .* sin(t) ...
    fp1 .* sin(t) + f .* cos(t) ];

rp2 = [ ( fp2 - f ) .* cos(t) - 2 * fp1 .* sin(t) ...
        ( fp2 - f ) .* sin(t) + 2 * fp1 .* cos(t) ];
    
% compute the Jacobian

Jacobian = sqrt( rp1(:,1).^2 + rp1(:,2).^2 );
     
% compute the curvature

kappa = ( rp1(:,1) .* rp2(:,2) - rp1(:,2) .* rp2(:,1) ) ...
    ./ Jacobian.^3;

% compute the unit tangent vector

tangent = [ rp1(:,1)./Jacobian rp1(:,2)./Jacobian ];

% compute the unit normal vector

normal = [ tangent(:,2) -tangent(:,1) ];   

end %return;