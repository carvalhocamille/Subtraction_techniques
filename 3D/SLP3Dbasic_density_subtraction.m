%% SLP3Dbasic_density_subtraction.m
%
% Computes the close evaluation of the exterior Neumann problem in 3D for 
% Laplace's equation using boundary integral methods to find the density of
% the single layer potential (SLP). 

% The SLP is computed for points close to the boundary on a body-fitted
% grid using 2 different methods.
%
% METHOD 1: a three-step method, in rotated coordinates that fixes the
% nearly singular point at the north pole (see Section 4 of Carvalho,
% Khatri, and Kim (2018) for details of the method implemented here.)
%
% METHOD 2: three-step method + density subtraction with linear function (DSL)
%
% This code calls the functions: 
%   For the thre-step method: Courtesy of A. D. Kim   
%       GaussLegendre.m: defines the wieghts and abscissas of Gauss-Legendre
%       quadrature
%       SphereRotation.m: rotates the coordinate system
%       ComputeSurface.:m: defines the surface's features
%       ComputeSphericalHarmonic.m: computes the spherical harmonics at the
%       grid points
%       Compute_Ymn.m: computes the spherical harmonics of degree n and order m
%   ComputeHarmonicFunction.m: computes the exact solution
%   ComputeSLPrho.m: computes the density (solve the BIE associated to the problem)
%
% Details of the method can be found in the manuscript 
% 'Layer potentials identities and subtraction techniques' Carvalho (2020).
% Written by C. Carvalho on 4/2/2020. 
% License MIT (see Readme.md)

clear all;
close all;
%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the order of the quadrature rules

N = 16;
M = 2*N;

%% compute the Gauss-Legendre quadrature rule in s

[ s, ws ] = GaussLegendre( N );
s  = 0.5 * pi * ( s + 1 );
ws = 0.5 * pi * ws;

%% compute the periodic trapezoid rule in t

dt = 2 * pi / M;
t  = -pi : dt : pi - dt;

%% compute mesh grid

[ S, T ] = meshgrid( s, t );

% stretch S and T into long column vectors

Svec = S(:);
Tvec = T(:);

%% set the order of the spherical harmonics approximation

n_order = N;

%% compute the DLP density spherical harmonics expansion coefficients

rho_nm = ComputeSLPrho( n_order );

disp( 'Computed density' );

%% plot the density coefficients to check resolution

figure; 
semilogy( (1:n_order^2), abs(rho_nm) )
ylabel( '$|\mu_{nm}|$', 'Interpreter', 'LaTeX' );


%% set the values of epsilon 

epsilon = [ 1e-8 1e-7 1e-6 1e-5 1e-4 ( 0.001 : 0.005 : 1 ) ];   %ellipse2, ellipse3-xz, mushroom
epsilon = fliplr(epsilon);

%% allocate memory for the evaluation of the asymptotic approximation

x1      = zeros( M * N, length(epsilon) );
x2      = zeros( M * N, length(epsilon) );
x3      = zeros( M * N, length(epsilon) );
approx1 = zeros( M * N, length(epsilon) );
approx2 = zeros( M * N, length(epsilon) );
exact   = zeros( M * N, length(epsilon) );
error1  = zeros( M * N, length(epsilon) );
error2  = zeros( M * N, length(epsilon) );

%% loop over all ystar values

indx = 8; % fix the value to just one ystar set by indx

for j = indx : indx %1 : 1 % N * M

    %% set the value of theta0 and phi0
    
    theta0 = Svec(j);
    phi0   = Tvec(j);

    %% compute the boundary data at (theta0,phi0)
    
    [ ~, ~, ystar, nustar, ~ ] = ComputeSurface( 0, 0, theta0, phi0 );

    %% compute the computational grid
    
    [ theta, phi, y, nu, J ] = ComputeSurface( theta0, phi0, Svec, Tvec );
    
    %% compute the spherical harmonics matrix

    Ynm     = ComputeSphericalHarmonics( n_order, theta, phi );
    Ynmstar = ComputeSphericalHarmonics( n_order, theta0, phi0 );

    %% compute the density

    rho     = Ynm * rho_nm;
    rhostar = Ynmstar * rho_nm;
    
    %% loop over all epsilon values
    
    for i = 1 : length( epsilon )

        %% compute the evaluation point

        x1(j,i) = ystar(1) + epsilon(i) * nustar(1);
        x2(j,i) = ystar(2) + epsilon(i) * nustar(2);
        x3(j,i) = ystar(3) + epsilon(i) * nustar(3);
                
        %% compute the exact solution
        
        exact(j,i) = ComputeHarmonicFunction( x1(j,i), x2(j,i), x3(j,i) );
        
        % compute the Y-vector
        
        Y1 = x1(j,i) - y(:,1);
        Y2 = x2(j,i) - y(:,2);
        Y3 = x3(j,i) - y(:,3);
        
        YX1 = y(:,1)-ystar(1);
        YX2 = y(:,2)-ystar(2);
        YX3 = y(:,3)-ystar(3);
        
        YXtot  = sqrt(YX1.^2 + YX2.^2 + YX3.^2);
        indexstar  = find(YXtot == min(YXtot) ); 
        
        % compute the length of Y
        
        Ylength = sqrt( Y1.^2 + Y2.^2 + Y3.^2 );
            
        % compute the kernel
        
        K = 0.5 * J .* sin(Svec) ./ Ylength;
        
        Kdlp = 0.5 * J .* sin(Svec) .* ...
            ( ( nu(:,1) .* Y1  + nu(:,2) .* Y2  + nu(:,3) .* Y3 ) ...
            ./ Ylength.^3 );
     
        % METHOD 1: three-step method
        
        Fbar1 = sum( reshape( K .* rho, M, N ), 1 ) / M;
        
        %METHOD 2: Density subtraction technique with lienar function
        
        usol1 = nustar(1) .* y(:,1) + nustar(2) .* y(:,2) + nustar(3) .* y(:,3);
        usol1star = nustar(1) .* ystar(1) + nustar(2) .* ystar(2) + nustar(3) .* ystar(3);
        dusol1 = nustar(1) .* nu(:,1)  + nustar(2) .* nu(:,2)  + nustar(3) .* nu(:,3);
        
        Fbar2 = sum( reshape( K .* rho .* (1 - dusol1) , M, N ), 1 ) / M ...
              + sum( reshape( K .* (rho - rhostar) .* dusol1 , M, N ), 1 ) / M ...
              + rhostar * sum( reshape( Kdlp .* (usol1 - usol1star) , M, N ), 1 ) / M ;
        
        
        %% compute the single-layer potential

        approx1(j,i) = real( sum( Fbar1 .* ws ) );
        approx2(j,i) = real( sum( Fbar2 .* ws ) );

        %% compute the error
        
        error1(j,i) = abs( exact(j,i) - approx1(j,i) );
        error2(j,i) = abs( exact(j,i) - approx2(j,i) );
       
    end

    %% output progress to the user

    disp( [ '   j = ', num2str(j), ' out of ', num2str(N*M) ] );
    
end

%% plot the error at a fixed point on the sphere

figure;
semilogx( epsilon, exact(indx,:), epsilon, approx1(indx,:), '--', ...
    epsilon, approx2(indx,:), '-.' );
xlabel( '$\epsilon$', 'Interpreter', 'LaTeX' );
ylabel( '$u( y^{\star} - \epsilon \nu^{\star} )$', 'Interpreter', 'LaTeX' );
legend( 'exact', 'method', 'DSL', 'Location', 'southwest' );

figure;
loglog( epsilon, error1(indx,:), 'o-', ...
    epsilon, error2(indx,:), 's-');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'error at point B', 'Interpreter', 'LaTeX' );
ylim([1e-16 1])
grid
legend( 'Method', 'DSL');
polyfit( log(epsilon), log(error1(indx,:)), 1 )
polyfit( log(epsilon), log(error2(indx,:)), 1 )

