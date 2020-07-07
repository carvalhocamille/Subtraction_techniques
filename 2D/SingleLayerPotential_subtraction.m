%%SingleLayerPotential.m
%
% This code computes the solution of the exterior Neumann problem for 
% Laplace's equation using boundary integral methods to find the density of
% the single layer potential (SLP). 
%
% The boundary integral equation is solved using the periodic trapezoid
% rule (PTR).
%
% The SLP is computed for points close to the boundary on a body-fitted
% grid using three different methods.
%
% METHOD 1: Native PTR that uses the same PTR used to solve the boundary
%           integral equation.
%
% METHOD 2: Density subtraction using linear function (DSL) and PTR
%
% METHOD 2: Density subtraction using Green-based function (DSG) and PTR
%
% This code calls the function, BoundaryCurve.m, to define the boundary of
% the domain.
%
% Details of the method can be found in the manuscript 
% 'Layer potentials identities and subtraction techniques' Carvalho (2020). 
% By C. Carvalho 07/2020
% License MIT (see Readme.md)

clear;
close all;

%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% compute the t grid

M  = 256;
dt = 2.0 * pi / M;
t  = ( 0 : dt : 2 * pi - dt ).';
kt = fftshift( - M / 2 : M / 2 - 1 ).';

%% compute the boundary curve and its properties

[ rB, tangent, normal, Jacobian, kappa ] = BoundaryCurve( t );

%% compute the Neumann data associated with the harmonic function x / r.

xbdy     = rB(:,1);
x0       = 0.1;
ybdy     = rB(:,2);
y0       = 0.4;
rsquared = ( ( xbdy - x0 ).^2 + ( ybdy - y0 ).^2 ).^2;
grad_x   = ( -( xbdy - x0 ).^2 + ( ybdy - y0 ).^2 ) ./ rsquared;
grad_y   = ( - 2 * ( xbdy - x0 ) .* ( ybdy - y0 ) ) ./ rsquared;
g        = normal(:,1) .* grad_x + normal(:,2) .* grad_y;

%% solve the boundary integral equation using the PTR

K = zeros( M );

for i = 1 : M
    
    K(i,i) = - 0.25 * kappa(i) / pi * Jacobian(i);
    
    for j = 1 : M
        
        if i ~= j
    
            K(i,j) = normal(i,1) * ( rB(i,1) - rB(j,1) ) ...
                + normal(i,2) * ( rB(i,2) - rB(j,2) );
        
            K(i,j) = - 0.5 / pi * K(i,j) * Jacobian(j) ...
                / ( ( rB(i,1) - rB(j,1) )^2 + ( rB(i,2) - rB(j,2) )^2 );
            
        end
        
    end
    
end

mu = ( K * dt - 0.5 * eye(M) ) \ g;

%% test solvability condition: 
%  the integral of the density over the boundary must be zero identically

disp( ' ' );
disp( ['   Solvability condition test (should be zero): ', ...
    num2str( sum( mu .* Jacobian ) * dt ) ] );
disp( ' ' );

%% compute the minimum and maximum curvature

kappa_max = max( abs(kappa) );

%% compute the epsilon grid points

ngrid = 200;
rgrid = linspace( 0, 1 / kappa_max, ngrid + 1 );
rgrid = rgrid(2:ngrid+1);

% rgrid = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 ( 0.001 : 0.005 : 1 ) ];  
% rgrid = fliplr(rgrid);
% ngrid  = length(rgrid);
%% allocate memory for the solutions

xplot = zeros( M+1, ngrid );
yplot = zeros( M+1, ngrid );
uplot = zeros( M+1, ngrid );
v1plot = zeros( M+1, ngrid );
v2plot = zeros( M+1, ngrid );
wplot = zeros( M+1, ngrid );

%% evaluate the solutions

for it = 1 : M
    
    for ir = 1 : ngrid
        
        %% compute the target point
        
        xplot(it,ir) = rB(it,1) + rgrid(ir) * normal(it,1);
        yplot(it,ir) = rB(it,2) + rgrid(ir) * normal(it,2);
        
        %% compute the native PTR kernel for the SLP
        
        ydiff  = [ xplot(it,ir) - rB(:,1) yplot(it,ir) - rB(:,2) ];
        distance  = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
        kernel = - 0.5 / pi .* log( distance );
                
        costheta  = sum( normal .* ydiff, 2 ) ./ distance;         
        kernel_dlp = 0.5 / pi * costheta ./ distance;

        %% METHOD 1 : native PTR 
        
        uplot(it,ir) = sum( kernel .* mu .* Jacobian ) * dt;

  
        %% METHOD 2: DSL
        nynx = normal(it,1) .* normal(:,1) + normal(it,2) .* normal(:,2);
        usol1 = normal(it,1) .* rB(:,1) + normal(it,2) .* rB(:,2);
        dusol1 = nynx;
        
        v1plot(it,ir) = sum( kernel .* (1- dusol1) .* mu .* Jacobian ) * dt  ...
            +  sum( kernel .* (mu - mu(it)) .* dusol1 .* Jacobian ) * dt  ...
            + mu(it) *  sum( kernel_dlp .* (usol1 - usol1(it)) .* Jacobian ) * dt ;
        
        %% METHOD 3: DSG
        yx  = [ rB(:,1)-rB(it,1)-normal(it,1), rB(:,2)-rB(it,2)-normal(it,2)];
        distanceyx  = sqrt( yx(:,1).^2 + yx(:,2).^2 );
        costhetayx  = sum( normal .* yx, 2 ) ./ distanceyx;
        
        usol2 = - log( distanceyx );
        dusol2 = - costhetayx ./ distanceyx;
                        
        v2plot(it,ir) = sum( kernel .* ( 1- dusol2) .* mu .* Jacobian ) * dt  ...
                    +  sum( kernel .* (mu - mu(it)) .* dusol2 .* Jacobian ) * dt  ...
                    + mu(it) *  sum( kernel_dlp .* (usol2 - usol2(it)) .* Jacobian ) * dt ;
        %% evaluate the exact solution
        
        wplot(it,ir) = ( xplot(it,ir) - x0 ) ...
            ./ ( ( xplot(it,ir) - x0 )^2 + ( yplot(it,ir) - y0 )^2 ); 
        
    end
       
end

%% report the maximum errors

disp( ['   Maximum absolute error made by native PTR: ', ...
    num2str( max( max( abs( uplot - wplot ) ) ) ) ] );

disp( ' ' );

disp( ['   Maximum absolute error made by DSL: ', ...
    num2str( max( max( abs( v1plot - wplot ) ) ) ) ] );

disp( ' ' );

disp( ['   Maximum absolute error made by DSG: ', ...
    num2str( max( max( abs( v2plot - wplot ) ) ) ) ] );

disp( ' ' );

%% impose periodicity

mu(M+1)      = mu(1);
xplot(M+1,:) = xplot(1,:);
yplot(M+1,:) = yplot(1,:);
uplot(M+1,:) = uplot(1,:);
v1plot(M+1,:) = v1plot(1,:);
v2plot(M+1,:) = v2plot(1,:);
wplot(M+1,:) = wplot(1,:);

%% Plot of the errors
% METHOD 1
figure; 
surf(xplot, yplot,log10( abs( uplot - wplot ) ) );
shading interp; view(2); colorbar; caxis([-16 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2  ,1),rB(M/2,2),'xk','LineWidth',2);
title( 'Error PTR');
axis tight;
axis square;
% METHOD 2
figure;
surf(xplot, yplot,log10( abs( v1plot - wplot ) ) );
shading interp; view(2); colorbar; caxis([-16 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2 +1  ,1),rB(M/2 +1 ,2),'xk','LineWidth',2);
title( 'Error DSL');
axis tight;
axis square;
% METHOD 3
figure;
surf(xplot, yplot,log10( abs( v2plot - wplot ) ) );
shading interp; view(2); colorbar; caxis([-16 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2 +1  ,1),rB(M/2 +1 ,2),'xk','LineWidth',2);
title( 'Error DSG');
axis tight;
axis square;
%% Log-log plot along the noraml
figure;
it0 = M/2 +1;
loglog( (rgrid), (abs( uplot(it0,:) - wplot(it0,:)) ), 'r', ...
      (rgrid), (abs( v1plot(it0,:) - wplot(it0,:) )), 'b--x',...
      (rgrid), (abs( v2plot(it0,:) - wplot(it0,:) )), 'm--.');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point A ' );
ylim([1.e-15 1])
grid on;
legend( 'PTR','DSL','DSG');

figure;
it0 = M/2+30;
loglog( (rgrid), (abs( uplot(it0,:) - wplot(it0,:)) ), 'r', ...
      (rgrid), (abs( v1plot(it0,:) - wplot(it0,:) )), 'b--x',...
      (rgrid), (abs( v2plot(it0,:) - wplot(it0,:) )), 'm--.');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point B ' );
grid on;
ylim([1.e-15 1])
legend( 'PTR','DSL','DSG');




