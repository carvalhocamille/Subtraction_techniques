%% SLP3Dbasic_density_subtraction.m
%
% Computes the close evaluation of the exterior Neumann problem in 3D for 
% Laplace's equation using boundary integral methods to find the density of
% the single layer potential (SLP). 

% The SLP is computed for points close to the boundary on a body-fitted
% grid using different representations.
%
% We use a three-step method, in rotated coordinates that fixes the
% nearly singular point at the north pole (see Section 4 of Carvalho,
% Khatri, and Kim (2018) for details of the method implemented here.)
%
% METHOD 0: standard representation with the three-step method
%
% METHOD 1: Linear function density subtraction with the three-step method
%
% METHOD 2: Green's function density subtraction with the three-step method
%
% METHOD 3: Quadratic difference density subtraction with the three-step method
%
% METHOD 4: Quadratic product density subtraction with the three-step method
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
% 'Modified representations for the close evaluation problem' Carvalho (2021). 
% By C. Carvalho 09/2021 
% License MIT (see Readme.md)

clear all;
close all;
%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the order of the quadrature rules
NN = [8 16];
epsilon = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];

errorN0 = zeros(length(NN),length(epsilon));
errorN1 = zeros(length(NN),length(epsilon));
errorN2 = zeros(length(NN),length(epsilon));
errorN3 = zeros(length(NN),length(epsilon));
errorN4 = zeros(length(NN),length(epsilon));
for ss = 1: length(NN)
    N = NN(ss);
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
    
    %figure;
    %semilogy( (1:n_order^2), abs(rho_nm) )
    %ylabel( '$|\mu_{nm}|$', 'Interpreter', 'LaTeX' );
    
    
    %% set the values of epsilon
    
    epsilon = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]; %[ 1e-8 1e-7 1e-6 1e-5 1e-4 ( 0.001 : 0.005 : 1 ) ];   %ellipse2, ellipse3-xz, mushroom
    epsilon = fliplr(epsilon);
    
    %% allocate memory for the evaluation of the asymptotic approximation
    
    x1      = zeros( M * N, length(epsilon) );
    x2      = zeros( M * N, length(epsilon) );
    x3      = zeros( M * N, length(epsilon) );
    approx0 = zeros( M * N, length(epsilon) );
    approx1 = zeros( M * N, length(epsilon) );
    approx2 = zeros( M * N, length(epsilon) );
    approx3 = zeros( M * N, length(epsilon) );
    approx4 = zeros( M * N, length(epsilon) );
    exact   = zeros( M * N, length(epsilon) );
    error0  = zeros( M * N, length(epsilon) );
    error1  = zeros( M * N, length(epsilon) );
    error2  = zeros( M * N, length(epsilon) );
    error3  = zeros( M * N, length(epsilon) );
    error4  = zeros( M * N, length(epsilon) );
    
%% loop over all ystar values
    indx = [floor(N*M/16), floor(N*M/8)];
    pts_names = ['A','B'];  
    for lll = 1: 2 %1 : 1 % N * M
        j = indx(lll);
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
            
            %% METHOD 0: standard representation            
           Fbar0 = sum( reshape( K .* rho, M, N ), 1 ) / M;
            
            %% METHOD 1: linear density subtraction            
            usol1 = nustar(1) .* y(:,1) + nustar(2) .* y(:,2) + nustar(3) .* y(:,3);
            usol1star = nustar(1) .* ystar(1) + nustar(2) .* ystar(2) + nustar(3) .* ystar(3);
            dusol1 = nustar(1) .* nu(:,1)  + nustar(2) .* nu(:,2)  + nustar(3) .* nu(:,3);
            
            Fbar1 = sum( reshape( K .* rho .* (1 - dusol1) , M, N ), 1 ) / M ...
                + sum( reshape( K .* (rho - rhostar) .* dusol1 , M, N ), 1 ) / M ...
                + rhostar * sum( reshape( Kdlp .* (usol1 - usol1star) , M, N ), 1 ) / M ;
            
            %% METHOD 2: Green's function density subtraction
            YXstar = sqrt((y(:,1)-ystar(1)-nustar(1)).^2 + (y(:,2)-ystar(2)-nustar(2)).^2 + (y(:,3)-ystar(3)-nustar(3)).^2);
            
            usol2 = 1 ./ YXstar ;
            usol2star = 1 ./ sqrt(nustar(1).^2 + nustar(2).^2 + nustar(3).^2) ;
            dusol2 = -( ( nu(:,1) .* (y(:,1)-ystar(1)-nustar(1)) + nu(:,2) .* (y(:,2)-ystar(2)-nustar(2))  + nu(:,3) .* (y(:,3)-ystar(3)-nustar(3)))./ YXstar.^3 );
            
            Fbar2 = sum( reshape( K .* rho .* (1 - dusol2) , M, N ), 1 ) / M ...
                + sum( reshape( K .* (rho - rhostar) .* dusol2 , M, N ), 1 ) / M ...
                + rhostar * sum( reshape( Kdlp .* (usol2 - usol2star) , M, N ), 1 ) / M ;
            
            %% METHOD 3: quadratic difference density subtraction
            coef3 = nustar(1) .* ystar(1) - nustar(2) .* ystar(2) ;
            usol3 = ( (y(:,1).* y(:,1)) - (y(:,2).* y(:,2))) / (2 * coef3);
            usol3star = ( (ystar(1).* ystar(1)) - (ystar(2).* ystar(2))) / (2* coef3);
            dusol3 = (nu(:,1) .* y(:,1) - nu(:,2) .* y(:,2))/coef3;
            
            Fbar3 = sum( reshape( K .* rho .* (1 - dusol3) , M, N ), 1 ) / M ...
                + sum( reshape( K .* (rho - rhostar) .* dusol3 , M, N ), 1 ) / M ...
                + rhostar * sum( reshape( Kdlp .* (usol3 - usol3star) , M, N ), 1 ) / M ;
            
            %% METHOD 4: quadratic product density subtraction
            coef4 = nustar(1) .* (ystar(2) - 5) + nustar(2) .* (ystar(1) - 5) ;
            usol4 = (y(:,1)-5) .* (y(:,2)-5) / coef4;
            usol4star = (ystar(1) -5) .* (ystar(2)-5) / coef4;
            dusol4 = (nu(:,1) .* (y(:,2)-5) + nu(:,2) .* (y(:,1)- 5))/coef4;
            
            Fbar4 = sum( reshape( K .* rho .* (1 - dusol4) , M, N ), 1 ) / M ...
                + sum( reshape( K .* (rho - rhostar) .* dusol4 , M, N ), 1 ) / M ...
                + rhostar * sum( reshape( Kdlp .* (usol4 - usol4star) , M, N ), 1 ) / M ;
            
            %% compute the single-layer potential
            approx0(j,i) = real( sum( Fbar0 .* ws ) );
            approx1(j,i) = real( sum( Fbar1 .* ws ) );
            approx2(j,i) = real( sum( Fbar2 .* ws ) );
            approx3(j,i) = real( sum( Fbar3 .* ws ) );
            approx4(j,i) = real( sum( Fbar4 .* ws ) );
            
            %% compute the error
            error0(j,i) = abs( exact(j,i) - approx0(j,i) );
            error1(j,i) = abs( exact(j,i) - approx1(j,i) );
            error2(j,i) = abs( exact(j,i) - approx2(j,i) );
            error3(j,i) = abs( exact(j,i) - approx3(j,i) );
            error4(j,i) = abs( exact(j,i) - approx4(j,i) );
                       
            errorN0(ss,i) =error0(j,i);
            errorN1(ss,i) = error1(j,i);
            errorN2(ss,i) = error2(j,i);
            errorN3(ss,i) = error3(j,i);
            errorN4(ss,i) = error4(j,i);
            
        end
        
        %% output progress to the user
        
        disp( [ '   j = ', num2str(j), ' out of ', num2str(N*M) ] );
        
        
        %% plot the error at a fixed point on the sphere        
        f = figure;
        loglog( epsilon, error0(j,:), 'k--', ...
            epsilon, error1(j,:), '--x', ...
            epsilon, error2(j,:), '--o',...
            epsilon, error3(j,:), '--.',...
            epsilon, error4(j,:), '--+');
        xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
        ylabel( 'Absolute error','Interpreter', 'LaTeX' );
        title(['Error at point ',num2str(pts_names(lll))],'Interpreter', 'LaTeX');
        legend( 'V0', 'V1','V2','V3','V4','Interpreter', 'LaTeX','Location', 'northeastoutside');
        xlim( [epsilon(end) epsilon(1)] )
        ylim( [ 1e-16 1e0 ] )
        grid on;
        f.Position = [10 10 550 400];
        %         %saveas(f,['SLP3D_Error_VS_lsmall_4methods_N',num2str(N),'_point',num2str(pts_names(lll)),'.png'])
        
    end
end

%%
for jj = 1:length(epsilon)
    f = figure;
        loglog( NN, errorN0(:,jj), 'k--', ...
            NN, errorN1(:,jj), '--x', ...
            NN, errorN2(:,jj), '--o',...
            NN, errorN3(:,jj), '--.',...
            NN, errorN4(:,jj), '--+');
        xlabel( '$N$', 'Interpreter', 'LaTeX' );
        ylabel( 'Absolute error','Interpreter', 'LaTeX' );
        title(['Error at point B, $\ell$ = ', num2str(epsilon(jj))],'Interpreter', 'LaTeX');
        legend( 'V0', 'V1','V2','V3','V4','Interpreter', 'LaTeX','Location', 'northeastoutside');
        xlim( [NN(1) NN(end)] )
        ylim( [ 1e-16 1e0 ] )
        grid on;
        f.Position = [10 10 550 400];
        %saveas(f,['SLP3D_Error_VS_M_4methods_l',num2str(epsilon(jj)),'_pointMs16.png'])
end