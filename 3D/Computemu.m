%% Computemu.m
%
% This code computes the coefficients of the spherical harmonics expansion
% for the density of the double- and single-layer potential on a closed surface using
% the discretized Galerkin method. 
% We test with and without density subtraction on the boundary
%
% Written by C. Carvalho on 4/2/2020.
% License MIT (see Readme.md)

function [c_coeffs1, c_coeffs2] = Computemu( n_order, kw )

%% QUADRATURE RULES

N = n_order;
M = 2 * N;

% Gauss-Legendre/trapezoid product rule for theta (z = cos(theta))

[ z, wz ] = GaussLegendre( N );

% transform z quadrature rule for s quadrature rule

s  = 0.5 * pi * ( z' + 1 );
ws = 0.5 * pi * wz';

% trapezoid rules for t

dt   =  pi / N;
t    = -pi : dt : pi - dt;

%% COMPOSITE ANGLE ARRAYS

PHI   = repmat( t', N, 1 );
THETA = repelem( acos( z' ), M );
W_GL  = repelem( wz', M ) * dt;

T     = repmat( t', N, 1 );
S     = repelem( s, M );

%% SPHERICAL HARMONICS MATRIX

[ Y_GL, ~ ] = ComputeSphericalHarmonics( n_order, THETA, PHI );

%% COMPUTE THE DIRICHLET BOUNDARY DATA

[ ~, ~, x, ~, ~ ] = ComputeSurface( 0, 0, THETA, PHI );
[ f, ~, ~, ~ ]    = ComputeFunction( x(:,1), x(:,2), x(:,3), kw );

%% COMPUTE SPHERICAL HARMONICS EXPANSION COEFFICIENTS

f_coeffs = Y_GL' * diag( W_GL ) * f;

%% ALLOCATE MEMORY FOR MATRICES

Kmtx1 = zeros( N * M, n_order^2 );
Kmtx2 = zeros( N * M, n_order^2 );
Kmtx = zeros( N * M, n_order^2 );
%% SET UP THE WAITBAR

wb = waitbar( 0, 'Computing the density. Please wait... ' );

%% COMPUTE THE FIRST INNER PRODUCT: BIE OPERATOR APPLIED TO SPHERICAL HARMONICS

for j = 1 : N * M
    
    % update the waitbar for the impatient user
    
    waitbar( j/(N*M), wb );
    
    % set the values of theta0 and phi0
    
    theta0 = THETA(j);
    phi0   = PHI(j);
   
    % compute ystar
    
    [ ~, ~, ystar, nustar, ~ ] = ComputeSurface( 0, 0, theta0, phi0 );

    % compute the surface
    
    [ varTHETA, varPHI, y, nu, J ] = ComputeSurface( theta0, phi0, S, T );
    
    % compute the kernel for the BIE operator
    
    yd = [ ystar(1) - y(:,1), ystar(2) - y(:,2), ystar(3) - y(:,3) ];
    
    ydist = sqrt( yd(:,1).^2 + yd(:,2).^2 + yd(:,3).^2 );
    
    nu_nustar = ( nustar(1) .* nu(:,1)  + nustar(2) .* nu(:,2)  + nustar(3) .* nu(:,3) ) ;
    
    nu_x_y = ( nu(:,1) .* yd(:,1)  + nu(:,2) .* yd(:,2)  + nu(:,3) .* yd(:,3) ) ;
    
    nustar_x_y = ( nustar(1) .* yd(:,1)  + nustar(2) .* yd(:,2)  + nustar(3) .* yd(:,3) ) ;
    
    
    SLP =  0.5 * J .* exp(1i * kw * ydist ) ./ ydist ;
    
    DLP =  (1./ydist - 1i*kw).* nu_x_y./ ydist  .* SLP;
    
    kernel = (DLP - 1i*kw * SLP);
    
    Kpw =  DLP - 1i*kw .* nu_nustar .* SLP;
    
    %compute plane wave
    PW = exp(- 1i * kw * nustar_x_y );
    


    % compute the integral operation
    
    k = 0;

    for n = 0 : n_order - 1
        
        for m = -n : n
            
            k = k + 1;
            
            % compute the function to be integrated
            
            % ktemp1 = ( Compute_Ynm( n, m, varTHETA, varPHI ) ...
            %     - Compute_Ynm( n, m, theta0, phi0 ) ) .* kernel .* sin(S);
                                    
            ktemp01 =  ( Compute_Ynm( n, m, varTHETA, varPHI )...
                - Compute_Ynm( n, m, theta0, phi0 ).*PW ).* (Kpw).* sin(S);
            
            ktemp02 = 1i * kw * (1 - nu_nustar) .* Compute_Ynm( n, m, varTHETA, varPHI).* SLP.* sin(S) ;
            
            ktemp1 =  Compute_Ynm( n, m, varTHETA, varPHI ) .* kernel .* sin(S);
            
            % compute the integral in t
            
            ktemp11 = sum( reshape( ktemp01, M, N ), 1 ) / M;
            ktemp12 = sum( reshape( ktemp02, M, N ), 1 ) / M;
            
            ktemp2 = sum( reshape( ktemp1, M, N ), 1 ) / M;          
            % compute the integral in s
            
            Kmtx1(j,k) = ktemp11 * ws;
            Kmtx2(j,k) = ktemp12 * ws;
            Kmtx(j,k)  = ktemp2 * ws;
        end
        
    end  
    
end

% close the wait bar

close( wb );

%% COMPUTE THE SECOND INNER PRODUCT

K1 = Y_GL' * diag( W_GL ) * Kmtx1;
K2 = Y_GL' * diag( W_GL ) * Kmtx2;

K = Y_GL' * diag( W_GL ) * Kmtx;

%% SOLVE THE GALERKIN SYSTEM OF EQUATIONS

%c_coeffs = ( 0.5 * eye( n_order^2 ) + K ) \ f_coeffs;

c_coeffs1 = ( K1 - K2 ) \ f_coeffs;

c_coeffs2 = ( 0.5 * eye( n_order^2 ) + K ) \ f_coeffs;
return
