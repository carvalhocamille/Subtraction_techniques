%% helmholtz3D.m
%
% Computes the close evaluation of the double- and single- layer potential in 3D for the sound-soft 
% acoustic scattering problem using a density subtraction technique (SLP-DLP).
%
%The SLP-DLP is computed for points close to the boundary on a body-fitted
% grid using 2 different representations. We also test for 2 density computations:
% We use a three-step method, in rotated coordinates that fixes the
% nearly singular point at the north pole (see Section 4 of Carvalho,
% Khatri, and Kim (2018) for details of the method implemented here.)
%
% Galerkin0: Product Gaussian Quadrature (PGQ) method to solve the standard BIE
% Galerkin1 : PGQ method to solve the modified BIE
%
% METHOD 0: standard representation with a three step-method and Galerkin0 density 
%
% METHOD 1: plane wave density subtraction with a three step-method and Galerkin0 density 
%
% METHOD 2: standard representation with a three step-method and Galerkin1 density 
%
% METHOD 3: plane wave density subtraction with a three step-method and Galerkin1 density 
%
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
% 'Modified representations for the close evaluation problem' Carvalho (2021). 
% By C. Carvalho 09/2021 
% License MIT (see Readme.md)

clear all;
close all;
%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% set the order of the quadrature rules

%% set the order of the quadrature rules
NN = [8 16];
epsilon = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];

errorN0 = zeros(2,length(NN),length(epsilon));
errorN1 = zeros(2,length(NN),length(epsilon));
errorN2 = zeros(2,length(NN),length(epsilon));
errorN3 = zeros(2,length(NN),length(epsilon));
KKw =  [5 20]; % [1, 2, 3, 4, 5, 10,20,30,40,50,60,70,80,90,100];
Error0 = zeros(length(NN),length(KKw));
Error1 =  zeros(length(NN),length(KKw));
Error2 =  zeros(length(NN),length(KKw));
Error3 =  zeros(length(NN),length(KKw));
for ss = 1:length(NN)
    N = NN(ss);
    M = 2*N;
    for idx = 1: length(KKw)
        kw = KKw(idx);
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
        
        %% allocate memory for the evaluation of the asymptotic approximation
        
        x1      = zeros( M * N, length(epsilon) );
        x2      = zeros( M * N, length(epsilon) );
        x3      = zeros( M * N, length(epsilon) );
        approx0 = zeros( M * N, length(epsilon) );
        approx1 = zeros( M * N, length(epsilon) );
        approx2 = zeros( M * N, length(epsilon) );
        approx3 = zeros( M * N, length(epsilon) );
        exact   = zeros( M * N, length(epsilon) );
        error0  = zeros( M * N, length(epsilon) );
        error1  = zeros( M * N, length(epsilon) );
        error2  = zeros( M * N, length(epsilon) );
        error3  = zeros( M * N, length(epsilon) );
        
        %% loop over all ystar values
        indx = [floor(N*M/16), floor(N*M/4)]; % fix the value to just one ystar set by indx
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
                
                
                %% METHOD 0: standard representation with three-step method and Galerkin0
                Fbar0 = sum( reshape( K .* mu, M, N ), 1 ) / M;        
                %% METHOD 1: modified representation with three-step method and Galerkin0
                PW = mustar * exp(- 1i * kw * nustar_x_y );               
                Fbar1 = sum( reshape( Kpw .* ( mu - PW ), M, N ), 1 ) / M; 
                Fbar11 = sum( reshape( SLP .* ( 1i*kw * (1 - nu_nustar) ).* mu, M, N ), 1 ) / M;
                
                %% METHOD 2: standard representation three-step method and Galerkin1                
                Fbar2 = sum( reshape( K .* mu2, M, N ), 1 ) / M;       
                %% METHOD 3: modified representation with three-step method and Galerkin1
                PW2 = mustar2 * exp(- 1i * kw * nustar_x_y );
                Fbar3 = sum( reshape( Kpw .* ( mu2 - PW2 ), M, N ), 1 ) / M; 
                Fbar33 = sum( reshape( SLP .* ( 1i*kw * (1 - nu_nustar) ).* mu2, M, N ), 1 ) / M;
                
                %% compute the double-layer potential
                approx0(j,i) = ( sum( Fbar0 .* ws ) );
                approx1(j,i) = ( sum( Fbar1 .* ws ) - sum( Fbar11 .* ws ) );
                approx2(j,i) = ( sum( Fbar2 .* ws ) );
                approx3(j,i) = ( sum( Fbar3 .* ws ) - sum( Fbar33 .* ws ) );
                %% compute the error
                error0(j,i) = abs( exact(j,i) - approx0(j,i) );
                error1(j,i) = abs( exact(j,i) - approx1(j,i) );
                error2(j,i) = abs( exact(j,i) - approx2(j,i) );
                error3(j,i) = abs( exact(j,i) - approx3(j,i) );
                
                errorN0(lll,ss,i) = error0(j,i);
                errorN1(lll,ss,i) = error1(j,i);
                errorN2(lll,ss,i) = error2(j,i);
                errorN3(lll,ss,i) = error3(j,i);
                
            end
            
            %% output progress to the user
            
            disp( [ '   j = ', num2str(j), ' out of ', num2str(N*M) ] );
            
            f = figure;
            loglog(epsilon, error0(j,:), 'k', ...
                epsilon, error1(j,:), '--x', ...
                epsilon, error2(j,:), '--o',...
                epsilon, error3(j,:), '--.');
            xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
            ylabel( 'Absolute error','Interpreter', 'LaTeX' );
            title(['Error at point ',num2str(pts_names(lll))],'Interpreter', 'LaTeX');
            legend( 'V0', 'V1','V2','V3','Interpreter', 'LaTeX','Location', 'northeastoutside');
            xlim( [epsilon(1) epsilon(end)] );
            ylim( [ 1e-16 1e0 ] );
            grid on;
            f.Position = [10 10 550 400];
        end
        Error0(ss,idx) = max(max(error0));
        Error1(ss,idx) = max(max(error1));
        Error2(ss,idx) = max(max(error2));
        Error3(ss,idx) = max(max(error3));
        %end
    end
    
    %% Error with respect to the wavenumber k
    f = figure;
    loglog(KKw,(Error0(ss,:)), 'k', KKw, (Error1(ss,:)), '--x');
    xlabel('$k$', 'Interpreter', 'LaTeX');
    ylim([1e-5 1e2]);
    legend( 'V0', 'V1','Interpreter', 'LaTeX','Location', 'northeastoutside');
    ylabel( 'Maximum error','Interpreter', 'LaTeX' );
    title(['Maximum error at point with respect to k'],'Interpreter', 'LaTeX');
    xlim( [ KKw(1) KKw(end)] );
    grid on;
    f.Position = [10 10 550 400];
    %saveas(f,['Hel3D_Error_inf_VS_k_2methods_N',num2str(N),'_pointA_rescaled.png']);
end
%% Error with respect to N
for lll = 1:2
    for jj = 1:length(epsilon)
        f = figure;
        loglog( NN, errorN0(lll,:,jj), 'k--', ...
            NN, errorN1(lll,:,jj), '--x', ...
            NN, errorN2(lll,:,jj), '--o',...
            NN, errorN3(lll,:,jj), '--.');
        xlabel( '$N$', 'Interpreter', 'LaTeX' );
        ylabel( 'Absolute error','Interpreter', 'LaTeX' );
        title(['Error at point ',num2str(pts_names(lll)),' $\ell$ = ', num2str(epsilon(jj))],'Interpreter', 'LaTeX');
        legend( 'V0', 'V1','V2','V3','Interpreter', 'LaTeX','Location', 'northeastoutside');
        xlim( [NN(1) NN(end)] )
        ylim( [ 1e-16 1e0 ] )
        grid on;
        f.Position = [10 10 550 400];
    end
end
