%% helmholtz3D.m
%
% Computes the close evaluation of the double- and single- layer potential in 3D for the sound-soft 
% acoustic scattering problem using a density subtraction technique (SLP-DLP).
%
% The SLP-DLP is computed for points close to the boundary on a body-fitted
% grid using 2 different methods.
%
% METHOD 1: a three-step method, in rotated coordinates that fixes the
% nearly singular point at the north pole (see Section 4 of Carvalho,
% Khatri, and Kim (2018) for details of the method implemented here.)
%
% METHOD 2: three-step method + plane wave density subtraction (PWS)
%
% We also test for 2 density computations:
%
% BIE-Method: solve the boundary integral equation as a Galerkin
% approximation
%
% BIE-PWS : solve the boundary integral equation as a Galerkin
% approximation + plane wave subtraction
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
%   ComputeFunction.m: computes the exact solution
%   Computemu.m: computes the density (solve the BIE associated to the problem)
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

Kw = 5 ;
Error3D = zeros(length(Kw),1);
ErrorPWS = zeros(length(Kw),1);
for idx = 1: length(Kw)
    
    %% set the wavenumber
    kw = Kw(idx);
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
    
    %% compute the density spherical harmonics expansion coefficients
    
    [mu_nm, mu_nm2]= Computemu( n_order, kw );

    disp( 'Computed density' );
    %% plot the density coefficients to check resolution
    
    figure(1);
    semilogy( (1:n_order^2), abs(mu_nm) )
    ylabel( '$|\mu_{nm}|$', 'Interpreter', 'LaTeX' );
    
    %% set the values of epsilon
    
    epsilon = [ 1e-6 1e-5 1e-4 ( 0.001 : 0.005 : 1 ) ];
    
    %% allocate memory for the evaluation of the asymptotic approximation
    
    x1      = zeros( M * N, length(epsilon) );
    x2      = zeros( M * N, length(epsilon) );
    x3      = zeros( M * N, length(epsilon) );
    approx1 = zeros( M * N, length(epsilon) );
    approx2 = zeros( M * N, length(epsilon) );
    approx3 = zeros( M * N, length(epsilon) );
    approx4 = zeros( M * N, length(epsilon) );
    exact   = zeros( M * N, length(epsilon) );
    error1  = zeros( M * N, length(epsilon) );
    error2  = zeros( M * N, length(epsilon) );
    error3  = zeros( M * N, length(epsilon) );
    error4  = zeros( M * N, length(epsilon) );
    
    %% loop over all ystar values
    
    indx = 200; % fix the value to just one ystar set by indx
    
    for j = indx: indx %1 : 1 % N * M
        
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
        
        mu     = Ynm * mu_nm;
        mustar = Ynmstar * mu_nm;
        
        mu2     = Ynm * mu_nm2;
        mustar2 = Ynmstar * mu_nm2;
        %% loop over all epsilon values
        
        for i = 1 : length( epsilon )
            
            %% compute the evaluation point
            
            x1(j,i) = ystar(1) + epsilon(i) * nustar(1);
            x2(j,i) = ystar(2) + epsilon(i) * nustar(2);
            x3(j,i) = ystar(3) + epsilon(i) * nustar(3);
            
            %% compute the exact solution
            
            [exact(j,i),~,~,~] = ComputeFunction( x1(j,i), x2(j,i), x3(j,i), kw );
            
            % compute the Y-vector
            
            Y1 = x1(j,i) - y(:,1);
            Y2 = x2(j,i) - y(:,2);
            Y3 = x3(j,i) - y(:,3);
            
            % compute the length of Y
            
            Ylength = sqrt( Y1.^2 + Y2.^2 + Y3.^2 );
            
            nu_nustar =  nustar(1) .* nu(:,1)  + nustar(2) .* nu(:,2)  + nustar(3) .* nu(:,3) ;
            
            nu_x_y =  nu(:,1) .* Y1  + nu(:,2) .* Y2  + nu(:,3) .* Y3  ;
            
            nustar_x_y =  nustar(1) .* Y1  + nustar(2) .* Y2  + nustar(3) .* Y3 ;
            
            % compute the kernel
            
            SLP =  0.5 * J .* exp(1i * kw * Ylength ) ./ Ylength .* sin(Svec);
            
            DLP =  (1./Ylength - 1i*kw).*( nu_x_y ./ Ylength ) .* SLP;
            
            K =  DLP - 1i*kw * SLP;
            
            Kpw =  DLP - 1i*kw * nu_nustar .* SLP;
            
            %compute plane wave
            
            PW = mustar * exp(- 1i * kw * nustar_x_y );
            PW2 = mustar2 * exp(- 1i * kw * nustar_x_y );
            
            % METHOD 1: three-step method
            
            Fbar1 = sum( reshape( K .* mu, M, N ), 1 ) / M;         % using BIE-Method
            Fbar4 = sum( reshape( K .* mu2, M, N ), 1 ) / M;        % using BIE-PWS
            
            % METHOD 2: PWS
            Fbar2 = sum( reshape( Kpw .* ( mu - PW ), M, N ), 1 ) / M; % using BIE-Method
            Fbar3 = sum( reshape( SLP .* ( 1i*kw * (1 - nu_nustar) ).* mu, M, N ), 1 ) / M;
            
           
            Fbar5 = sum( reshape( Kpw .* ( mu2 - PW2 ), M, N ), 1 ) / M; % using BIE-PWS
            Fbar6 = sum( reshape( SLP .* ( 1i*kw * (1 - nu_nustar) ).* mu2, M, N ), 1 ) / M;
            
            %% compute the double-layer potential
            
            approx1(j,i) = ( sum( Fbar1 .* ws ) );
            approx2(j,i) = ( sum( Fbar2 .* ws ) - sum( Fbar3 .* ws ) );
            
            approx3(j,i) = ( sum( Fbar4 .* ws ) );
            approx4(j,i) = ( sum( Fbar5 .* ws ) - sum( Fbar6 .* ws ) );
            %% compute the error
            
            error1(j,i) = abs( exact(j,i) - approx1(j,i) );
            error2(j,i) = abs( exact(j,i) - approx2(j,i) );
            
            error3(j,i) = abs( exact(j,i) - approx3(j,i) );
            error4(j,i) = abs( exact(j,i) - approx4(j,i) );
            
        end
        
        %% output progress to the user
        
        disp( [ '   j = ', num2str(j), ' out of ', num2str(N*M) ] );
        
    end
    Error3D(idx) = max(max(error1));
    ErrorPWS(idx) = max(max(error2));
end
%% Error with respect to the wavenumber k 
figure; 
semilogy(Kw,(Error3D), 'r', Kw, (ErrorPWS), '--b');
xlabel('k');
legend('Method','PWS');


%% Log-Log plot of the error along the normal Method and PWS
figure;
loglog( epsilon, error3(indx,:), 'o-', ...
    epsilon, error4(indx,:), 'x-');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'error', 'Interpreter', 'LaTeX' );
grid
legend('Method', 'PWS' );
ylim([1e-6 1])
    
 %% Log-Log plot of the error along the normal using BIE-Method and BIE-PWS
 figure;
 loglog( epsilon, error2(indx,:), 'r' , ...
     epsilon, error4(indx,:), 'b--');
 xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
 ylabel( 'error point A ', 'Interpreter', 'LaTeX' );
 grid
 legend( 'BIE Method', 'BIE-PWS');
 ylim([1e-6 1])

%% Comparison densities
figure;
semilogy(1:M*N, abs(mu-mu2))
xlabel('Quadrature point')
ylabel('comparion densities')

