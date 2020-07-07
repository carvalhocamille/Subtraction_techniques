%%Helmholtz_subtraction_techniques.m
%
% This code computes the solution of the sound-soft scattering problem 
% using boundary integral methods to find the density of
% the single- and double-layer potential (SLP-DLP). 
%
% The boundary integral equation is solved using Kress quadrature (for the density) and the periodic trapezoid
% rule (PTR, for the solution).
%
% The SLP-DLP is computed for points close to the boundary on a body-fitted
% grid using 3 different methods.
%
% METHOD 1: Native PTR.
%
% METHOD 2: PTR + plane wave density subtraction
%
% We also test for 2 density computations:
%
% KRESS: Kress quadrature to compute the solution mu of the BIE
%
% PTRds : Periodic Trapezoid Rule and density subtraction in the BIE
%
% This code calls the function, BoundaryCurve.m, to define the boundary of
% the domain.
%
% Details of the method can be found in the manuscript 
% 'Layer potentials identities and subtraction techniques' Carvalho (2020). 
% By C. Carvalho 07/2020
% License MIT (see Readme.md)

clear all;
close all;

%% set the figure parameters

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',1.0,...
      'defaultlinelinewidth',2.0,'defaultpatchlinewidth',1.0); 

%% compute the theta grid

M  = 256;
dt = 2.0 * pi / M;
t  = ( 0 : dt : 2 * pi - dt ).';
kt = fftshift( -M / 2 : M / 2 - 1 ).';

%% compute the boundary curve and its properties

[ rB, tangent, normal, Jacobian, kappa ] = BoundaryCurve( t );

%% set the wavenumber
Kw = 5;
ErrorPTR = zeros(length(Kw),1);
ErrorPWS = zeros(length(Kw),1);
for idx = 1: length(Kw)
    k = Kw(idx);
    
    %% compute the Dirichlet data
    
    xbdy   = rB(:,1);
    ybdy   = rB(:,2);
    x0     =  0.2;
    y0     = 0.8;
    dist   = sqrt( ( xbdy - x0 ).^2 + ( ybdy - y0 ).^2 );
    % alpha = pi/4;
    % f      = exp(1i * k * (xbdy * cos(alpha) + ybdy * sin(alpha)));
    f = (1i / 4) * besselh( 0, k * dist );
    
    %% KRESS: solve the boundary integral equation: Kress's Nystrom method
    
    % compute Euler's constant    
    C = 0.57721566490153286;
    Const =  1i * 0.25 - 0.5 / pi * (log( 0.5 * k) + C );
    
    % Kress quadrature (Barnett implementation)
    n  = M / 2;
    m  = 1 : n - 1;
    Rj = -2 * pi * ifft( [0 1./m 1/n 1./m(end:-1:1)] );
    R  = toeplitz( [ Rj(1) Rj(end:-1:2)], Rj );
    
    % solve the boundary integral equation
    
    Kbdy = nan( M, M );
    L1   = nan( M, M );
    L2   = nan( M, M );
    M1   = nan( M, M );
    M2   = nan( M, M );
    
    for i = 1 : M
        
        for j = 1 : M
            
            rdiff    = rB(i,:) - rB(j,:);
            distance = sqrt( rdiff(1).^2 + rdiff(2).^2 );
            
            M1(i,j) = - 0.5 / pi * Jacobian(j) * besselj( 0, k * distance );
            
            if i == j
                
                L1(i,j) = 0;
                
                L2(i,j) = 0.5 / pi * kappa(j) * Jacobian(j);
                
                M2(i,j) = Jacobian(j) * ( 0.5 * 1i - C / pi - 0.5 / pi ...
                    * log( 0.25 * k^2 * Jacobian(j)^2 ) );
                
            else
                
                costheta = normal(j,:) * rdiff' / distance;
                
                logterm = log( 4 * sin( 0.5 * ( i - j ) * dt )^2 );
                
                L1(i,j) = 0.5 * k / pi * Jacobian(j) * costheta ...
                    * besselj( 1, k * distance );
                
                L2(i,j) = -0.5 * 1i * k * Jacobian(j) * costheta ...
                    * besselh( 1, k * distance ) - L1(i,j) * logterm;
                
                M2(i,j) = 0.5 * 1i * Jacobian(j) * besselh( 0, k * distance ) ...
                    - M1(i,j) * logterm;
                
            end
            
            Kbdy(i,j) = R(i,j) * ( L1(i,j) + 1i * k * M1(i,j) ) ...
                + (pi / n ) * ( L2(i,j) + 1i * k * M2(i,j) );
            
        end
        
    end
    
    % solve for the density
    mu = ( eye(M) - Kbdy ) \ (2 * f);
    
    %% PTRds: solve the boundary integral equation: PTR + density subtraction
    
    K2 = zeros(M,M);
    Kw2 = zeros(M,M);
    PW = zeros(M,M);
    PW2 = zeros(M,1);
    for i = 1 : M
        
        for j = 1 : M
            if i ~= j
                
                rdiffb    = rB(i,:) - rB(j,:);
                dist2 = sqrt( rdiffb(1).^2 + rdiffb(2).^2 );
                costhetab = normal(j,:) * rdiffb' / dist2;
                
                %Helmholtz kernels
                kernel_dlpb    =  1i * 0.25 *  k * costhetab * besselh( 1, k * dist2 ) * Jacobian(j) * dt;
                kernel_slpb    =  1i * 0.25 * besselh( 0, k * dist2 ) * Jacobian(j) * dt;
                kernelb    =  kernel_dlpb - 1i * k * kernel_slpb;
                
                %plane wave technique prep
                costheta2 = (normal(j,1) * normal(i,1) + normal(j,2) * normal(i,2) );
                PW(i,j) = exp( -1i * k * ( normal(i,1) * rdiffb(1) + normal(i,2) * rdiffb(2) ) );
                
                %kernel in the plane wave direction
                kerneldb    =  kernel_dlpb - 1i * k * costheta2 * kernel_slpb;
                
                Kw2(i,j) = 1i * k * (1 - costheta2) * kernel_slpb ;
                K2(i,j)  = kerneldb ;
                PW(i,j) = K2(i,j) * PW(i,j);
            end
            
        end
        PW2(i) = sum(PW(i,:));
        
    end

    Kbdy2 = ( K2 - (diag(PW2)) - Kw2  ) ;  
    mu2 = ( Kbdy2 ) \ f;
    %% Comparison densities
    figure;
    plot (1:M, real(mu),'r', 1:M, imag(mu),'b',1:M, real(mu2),'--r', 1:M, imag(mu2),'--b' );
    legend('real mu', 'imag mu','real mu2', 'imag mu2');
    
    figure;
    semilogy(1:M, (abs(mu(1:M) -mu2(1:M))));
    xlabel('Quadrature point')
    ylabel('comparison densities')
    ylim([1e-7 1e-4])
    
    disp( ['   Maximum absolute error between mu and mu2: ', ...
        num2str( max(max(abs(mu - mu2) ) ) ) ] );
    
    disp( ' ' );
    
    %% compute the maximum curvature
    
    kappa_max = max( abs( kappa ) );
    
    %% compute the local grid points
    
    ngrid = 100;
    rgrid = linspace( 0, 1 / kappa_max, ngrid + 1 );
    rgrid = rgrid(2:ngrid);
    
%     rgrid = [1e-6 1e-5 1e-4 (0.001:0.005:1) 1.2];
%     ngrid = length(rgrid);
    %% allocate memory for the solution
    
    xplot = zeros( M+1, ngrid - 1 );
    yplot = zeros( M+1, ngrid - 1 );
    u1plot = zeros( M+1, ngrid - 1 );
    u2plot = zeros( M+1, ngrid - 1 );
    u3plot = zeros( M+1, ngrid - 1 );
    u4plot = zeros( M+1, ngrid - 1 );
    u5plot = zeros( M+1, ngrid - 1 );
    exactplot = zeros( M+1, ngrid - 1 );
    
    %% evaluate the solutions

    for it = 1 : M
        
        for ir = 1 : ngrid - 1
            
            %% compute the target point
            
            xplot(it,ir) = rB(it,1) + rgrid(ir) * normal(it,1);
            yplot(it,ir) = rB(it,2) + rgrid(ir) * normal(it,2);
            
            %% compute the kernels
            
            ydiff     = [ xplot(it,ir) - rB(:,1) yplot(it,ir) - rB(:,2) ];
            distance  = sqrt( ydiff(:,1).^2 + ydiff(:,2).^2 );
            costheta  = ( normal(:,1) .* ( xplot(it,ir) - rB(:,1) ) ...
                + normal(:,2) .* ( yplot(it,ir) - rB(:,2) ) ) ./ distance;
            
            %Helmholtz kernels
            kernel_dlp    =  1i * 0.25 *  k * costheta .* besselh( 1, k * distance ) ;
            kernel_slp    =  1i * 0.25 * besselh( 0, k * distance ) ;
            kernel    =  kernel_dlp - 1i * k * kernel_slp;
            
            %% METHOD 1: native PTR 
            
            u1plot(it,ir) = sum( kernel .* Jacobian .* mu ) * dt;  % using KRESS
            u4plot(it,ir) = sum( kernel .* Jacobian .* mu2 ) * dt; % using PTRds
            
            %% METHOD 2: PTR + plane wave density subtraction
            
            nx = [ xplot(it,ir) / (xplot(it,ir)^2 + yplot(it,ir)^2)^(1/2)     yplot(it,ir) / (xplot(it,ir)^2 + yplot(it,ir)^2)^(1/2)];
            nxny = normal(:,1)*nx(1) + normal(:,2)*nx(2) ;
            nyx = (rB(:,1)-rB(it,1))*nx(1) + (rB(:,2)-rB(it,2))*nx(2);
            usol = exp(1i *k * nyx);
            dusol =1i *k * nxny .* usol;
            
            u2plot(it,ir) = sum( (kernel_dlp -  kernel_slp .* dusol ) .* Jacobian .* (mu - mu(it) ) ) * dt ...
                            + sum( (dusol - 1i*k ).* kernel_slp .* Jacobian .* mu )  * dt ...
                            + mu(it) * sum( (1 - usol ).* kernel_dlp .* Jacobian)  * dt;                            % using KRESS
            
            u5plot(it,ir) = sum( (kernel_dlp -  kernel_slp .* dusol ) .* Jacobian .* (mu2 - mu2(it) ) ) * dt ...
                            + sum( (dusol - 1i*k ).* kernel_slp .* Jacobian .* mu2 )  * dt ...
                            + mu2(it) * sum( (1 - usol ).* kernel_dlp .* Jacobian)  * dt;                           % using PTRds
            
            %% Exact solution
            
            dist0 = sqrt( ( xplot(it,ir) - x0 ).^2 + ( yplot(it,ir) - y0 ).^2 );
            exactplot(it,ir) = (1i / 4) * besselh( 0, k * dist0 );
            
        end
        
    end
    
    %% report the maximum errors
    
    disp( ['   Maximum absolute error made by METHOD 1: ', ...
        num2str( max( max( abs( exactplot - u1plot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by METHOD 2: ', ...
        num2str( max( max( abs( exactplot - u2plot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by METHOD 1 and mu2: ', ...
        num2str( max( max( abs( exactplot - u4plot ) ) ) ) ] );
    
    disp( ' ' );
    
     disp( ['   Maximum absolute error made by METHOD 2 and mu2: ', ...
         num2str( max( max( abs( exactplot - u5plot ) ) ) ) ] );
     
    disp( ' ' );
    
    ErrorPTR(idx) =  max( max( abs( exactplot - u1plot ) ) );
    ErrorPWS(idx) =  max( max( abs( exactplot - u2plot ) ) );
end
%% Error with respect to the wavenumber k 
figure; 
semilogy(Kw,(ErrorPTR), 'r', Kw, (ErrorPWS), '--b');
xlabel('k');
ylim([1e-5 1e5]);
legend('PTR','PWS');

%% impose periodicity

mu(M+1)      = mu(1);
mu2(M+1)     = mu2(1);
xplot(M+1,:) = xplot(1,:);
yplot(M+1,:) = yplot(1,:);
u1plot(M+1,:) = u1plot(1,:);
u2plot(M+1,:) = u2plot(1,:);
u3plot(M+1,:) = u3plot(1,:);
u4plot(M+1,:) = u4plot(1,:);
u5plot(M+1,:) = u5plot(1,:);

exactplot(M+1,:) = exactplot(1,:);

%% plot the errors

figure; 
surf(xplot, yplot,log10( abs( u1plot - exactplot ) ) );
shading interp; view(2); colorbar; caxis([-14 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2+1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
title( 'error PTR');
axis tight;
axis square;

figure;
surf(xplot, yplot,log10( abs( u2plot - exactplot ) ) );
shading interp; view(2); colorbar; caxis([-14 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2 +1  ,1),rB(M/2 +1 ,2),'xk','LineWidth',2);
title( 'error PWS');
axis tight;
axis square;

figure;

surf(xplot, yplot,log10( abs( u4plot - exactplot ) ) );
shading interp; view(2); colorbar; caxis([-14 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2 +1  ,1),rB(M/2 +1 ,2),'xk','LineWidth',2);
title( 'error PTR mu2');
axis tight;
axis square;

figure;
surf(xplot, yplot,log10( abs( u5plot - exactplot ) ) );
shading interp; view(2); colorbar; caxis([-14 1]);
hold on;
scatter(rB(M / 2 +30 ,1),rB(M/2 +30,2),'xk','LineWidth',2);
hold on;
scatter(rB(M / 2 +1  ,1),rB(M/2 +1 ,2),'xk','LineWidth',2);
title( 'error PWS mu2');
axis tight;
axis square;

%% Log-Log plot Error PTR-PWS
figure;
it0 = M/2 +1;
loglog( (rgrid), (abs( u1plot(it0,:) - exactplot(it0,:)) ), 'r', ...
      (rgrid), (abs( u2plot(it0,:) - exactplot(it0,:) )), '--bx');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point A ' );
grid on;
ylim([1e-20 1e5]);
yticks([1e-20 1e-15 1e-10 1e-5 1 1e5]);
legend( 'PTR','PWS');

figure;
it0 = M/2+30;
loglog( (rgrid), (abs( u1plot(it0,:) - exactplot(it0,:)) ), 'r', ...
      (rgrid), (abs( u2plot(it0,:) - exactplot(it0,:) )), '--bx');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point B ' );
grid on;
ylim([1e-20 1e5]);
yticks([1e-20 1e-15 1e-10 1e-5 1 1e5]);
legend( 'PTR','PWS');


%% Comparison mu, mu2 on PWS

figure;
it0 = M/2 +1;
loglog( (rgrid), (abs( u2plot(it0,:) - exactplot(it0,:)) ), 'b', ...
      (rgrid), (abs( u5plot(it0,:) - exactplot(it0,:) )), '--rx');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point A ' );
grid on;
ylim([1e-16 1]);
legend( 'BIE-Method','BIE-PWS');

figure;
it0 = M/2+30;
loglog( (rgrid), (abs( u2plot(it0,:) - exactplot(it0,:)) ), 'b', ...
      (rgrid), (abs( u5plot(it0,:) - exactplot(it0,:) )), '--rx');
xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
ylabel( 'absolute error at point B ' );
grid on;
ylim([1e-16 1]);
legend( 'BIE-Method','BIE-PWS');


