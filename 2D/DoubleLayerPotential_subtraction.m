%%DoubleLayerPotential.m
%
% This code computes the solution of the interior Dirichlet problem for 
% Laplace's equation using boundary integral methods to find the density of
% the double layer potential (DLP). 
%
% The boundary integral equation is solved using the periodic trapezoid
% rule (PTR).
%
% The DLP is computed for points close to the boundary on a body-fitted
% grid using three different methods.
%
% METHOD 0: standard representation with PTR 
%
% METHOD 1: Classic constant density subtraction with PTR
%
% METHOD 2: Quadratic difference density subtraction with PTR
%
% METHOD 3: Quadratic product density subtraction with PTR
%
%
% This code calls the function, BoundaryCurve.m, to define the boundary of
% the domain.
%
% Details of the method can be found in the manuscript 
% 'Modified representations for the close evaluation problem' Carvalho (2021). 
% By C. Carvalho 09/2021
% License MIT (see Readme.md)
clear;
close all;

%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% compute the theta grid

M  = 128;
dt = 2.0 * pi / M;
t  = ( 0 : dt : 2 * pi - dt ).';
kt = fftshift( -M / 2 : M / 2 - 1 ).';

%% compute the boundary curve and its properties

[ rB, tangent, normal, Jacobian, kappa ] = BoundaryCurve( t );

%% compute the Dirichlet data

x0   = 2;
y0   = 1.8;
dist = sqrt( ( rB(:,1) - x0 ).^2 + ( rB(:,2) - y0 ).^2 );
f    = - 0.5 / pi * log( dist );

%% solve the boundary integral equation

Kbdy = zeros( M );

for i = 1 : M
    
    Kbdy(i,i) = -0.25 * kappa(i) / pi;
      
    for j = 1 : M
        
        if i ~= j

            ydiff     = [ rB(i,1) - rB(j,1) rB(i,2) - rB(j,2) ];
            distance  = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
            costheta  = sum( normal(j,:) .* ydiff, 2 ) ./ distance;                           
            Kbdy(i,j) = 0.5 / pi * costheta / distance;
            
        end
        
    end
    
end

mu = ( Kbdy * diag( Jacobian ) * dt - 0.5 * eye( M ) ) \ f;

%% compute the Fourier series coefficients of mu * Jacobian

Uhat = fft( mu .* Jacobian ) / M;

%% compute the maximum curvature

kappa_max = max( abs( kappa ) );

%% compute the local grid points

% ngrid = 99;
% rgrid = linspace( 0, 3.0 / kappa_max, ngrid + 1 ); 
% rgrid = rgrid(2:ngrid);

rgrid = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 ( 0.001 : 0.005 : 1 ) ];  
rgrid = fliplr(rgrid);
ngrid  = length(rgrid);
%% allocate memory for the solution

xplot = zeros( M+1, ngrid - 1 );
yplot = zeros( M+1, ngrid - 1 );
uplot = zeros( M+1, ngrid - 1 );
v0plot = zeros( M+1, ngrid - 1 );
v1plot = zeros( M+1, ngrid - 1 );
v2plot = zeros( M+1, ngrid - 1 );
wplot = zeros( M+1, ngrid - 1 );

%% evaluate the solutions

for it = 1 : M
    
    for ir = 1 : ngrid 
    
        %% compute the target point        
        xplot(it,ir) = rB(it,1) - rgrid(ir) * normal(it,1);
        yplot(it,ir) = rB(it,2) - rgrid(ir) * normal(it,2);
        
        %% compute the native PTR kernel for the DLP        
        ydiff     = [ xplot(it,ir) - rB(:,1) yplot(it,ir) - rB(:,2) ];
        distance  = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
        costheta  = sum( normal .* ydiff, 2 ) ./ distance;         
        kernel    = 0.5 / pi * costheta ./ distance;
        kernel_slp    =  - 0.5 / pi * log(distance);
        
        %% evaluate the native PTR DLP
        
        uplot(it,ir) = sum( kernel .* Jacobian .* mu ) * dt;
        
        %% Method 0: classic constant density subtraction        
        v0plot(it,ir) = sum( kernel .* Jacobian .* (mu - mu(it)) ) * dt - mu(it);
        
        %% Method 1: quadratic difference density subtraction
        usol1 = 1 + (rB(:,1)-rB(it,1)).*(rB(:,1)-rB(it,1)) - (rB(:,2)-rB(it,2)).*(rB(:,2)-rB(it,2)) ; 
        dusol1 = 2 * normal(:,1) .*(rB(:,1)-rB(it,1)) -2 * normal(:,2) .*(rB(:,2)-rB(it,2));
        v1plot(it,ir) = -mu(it) * (usol1(it) + dusol1(it)) ...
            + sum( ( kernel ) .* Jacobian .* mu .* (1 - usol1) ) * dt ...
            + sum( ( kernel ) .* Jacobian .* (mu - mu(it)) .* usol1 ) * dt ...
            + mu(it) * sum( ( kernel_slp ) .* Jacobian .* (dusol1 - dusol1(it)) ) * dt;
        
        %% Method 2: quadratic product density subtraction
        usol2 = 1 + (rB(:,1)-rB(it,1)).*(rB(:,2)-rB(it,2))  ; 
        dusol2 =  normal(:,1) .*(rB(:,2)-rB(it,2)) + normal(:,2) .*(rB(:,1)-rB(it,1));
        v2plot(it,ir) = -mu(it) * (usol2(it) + dusol2(it)) ...
            + sum( ( kernel ) .* Jacobian .* mu .* (1 - usol2) ) * dt ...
            + sum( ( kernel ) .* Jacobian .* (mu - mu(it)) .* usol2 ) * dt ...
            + mu(it) * sum( ( kernel_slp ) .* Jacobian .* (dusol2 - dusol2(it)) ) * dt;
             
        %% evaluate the exact solution        
        dist = sqrt( ( xplot(it,ir) - x0 ).^2 + ( yplot(it,ir) - y0 ).^2 );
        wplot(it,ir) = - 0.5 / pi * log( dist );
        
    end

end

%% report the maximum errors

disp( ['   Maximum absolute error made by native PTR: ', ...
    num2str( max( max( abs( uplot - wplot ) ) ) ) ] );

disp( ' ' );

disp( ['   Maximum absolute error made by constant density sub: ', ...
    num2str( max( max( abs( v0plot - wplot ) ) ) ) ] );

disp( ' ' );

disp( ['   Maximum absolute error made by quadratic difference density sub: ', ...
    num2str( max( max( abs( v1plot - wplot ) ) ) ) ] );

disp( ' ' );


disp( ['   Maximum absolute error made by quadratic product density sub: ', ...
    num2str( max( max( abs( v2plot - wplot ) ) ) ) ] );

disp( ' ' );


%% impose periodicity

mu(M+1)      = mu(1);
xplot(M+1,:) = xplot(1,:);
yplot(M+1,:) = yplot(1,:);
uplot(M+1,:) = uplot(1,:);
v0plot(M+1,:) = v0plot(1,:);
v1plot(M+1,:) = v1plot(1,:);
v2plot(M+1,:) = v2plot(1,:);
wplot(M+1,:) = wplot(1,:);

%% plot a slice of the errors on a semilog scale
points = [M/8, M/2+1, M/2 + 30]; 
pts_names = ['A','B','C'];

for idx = 1: 3
    it0 = points(idx);
    f = figure;
    loglog( rgrid(1:end), abs( uplot(it0,:) - wplot(it0,:) ), 'k--', ...
        rgrid(1:end), abs( v0plot(it0,:) - wplot(it0,:) ), '--x', ...
        rgrid(1:end), abs( v1plot(it0,:) - wplot(it0,:) ), '--.' ,...
        rgrid(1:end), abs( v2plot(it0,:) - wplot(it0,:) ), '--o');
    xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
    ylabel( 'absolute error' );
    legend( 'Standard', 'V0','V1', 'V2','Location', 'northeastoutside');
    title(['Error at point ',num2str(pts_names(idx))],'Interpreter','LaTeX');
    xlim( [ 1.e-10 1] )
    ylim( [ 1e-15 1e5 ] )
    grid on;
    f.Position = [10 10 550 400]; 
end
