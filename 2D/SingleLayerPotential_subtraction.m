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
% METHOD 0: standard representation with PTR 
%
% METHOD 1: Linear function density subtraction with PTR
%
% METHOD 2: Green's function density subtraction with PTR
%
% METHOD 3: Quadratic difference density subtraction with PTR
%
% METHOD 4: Quadratic product density subtraction with PTR
%
% METHOD 5: Exponential and Sin density subtraction with PTR
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

%% compute the t grid
MM = [128 256 512];
rgrid = [ 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
rgrid = fliplr(rgrid);
ngrid  = length(rgrid);
Err0 = zeros( length(MM), ngrid);
Err1 = zeros( length(MM), ngrid);
Err2 = zeros( length(MM), ngrid);
Err3 = zeros( length(MM), ngrid);
Err4 = zeros( length(MM), ngrid);
Err5 = zeros( length(MM), ngrid);
ErrDIM = zeros( length(MM), ngrid);
for idx = 1:length(MM)
    M  = MM(idx);
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
    Uhat = fft( mu .* Jacobian ) / M;
    %% test solvability condition:
    %  the integral of the density over the boundary must be zero identically
    disp( ' ' );
    disp( ['   Solvability condition test (should be zero): ', ...
        num2str( sum( mu .* Jacobian ) * dt ) ] );
    disp( ' ' );
    
    %% compute the minimum and maximum curvature
    
    kappa_max = max( abs(kappa) );
    
    %% allocate memory for the solutions
    
    xplot = zeros( M+1, ngrid );
    yplot = zeros( M+1, ngrid );
    v0plot = zeros( M+1, ngrid );
    v1plot = zeros( M+1, ngrid );
    v2plot = zeros( M+1, ngrid );
    v3plot = zeros( M+1, ngrid );
    v4plot = zeros( M+1, ngrid );
    v5plot = zeros( M+1, ngrid );
    DIMplot = zeros( M+1, ngrid );
    wplot = zeros( M+1, ngrid );
    %% evaluate the solutions
    it0 = floor(M/8);
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
            
            % METHOD 0 : standard representation
            v0plot(it,ir) = sum( kernel .* mu .* Jacobian ) * dt;
            
            
            %% METHOD 1: linear density subtraction
            nynx = normal(it,1) .* normal(:,1) + normal(it,2) .* normal(:,2);
            usol1 = normal(it,1) .* rB(:,1) + normal(it,2) .* rB(:,2);
            dusol1 = nynx;
            
            v1plot(it,ir) = sum( kernel .* (1- dusol1) .* mu .* Jacobian ) * dt  ...
                +  sum( kernel .* (mu - mu(it)) .* dusol1 .* Jacobian ) * dt  ...
                + mu(it) *  sum( kernel_dlp .* (usol1 - usol1(it)) .* Jacobian ) * dt ;
            
            %% METHOD 2: Green's function density subtraction
            yx  = [ rB(:,1)-rB(it,1)-normal(it,1), rB(:,2)-rB(it,2)-normal(it,2)];
            distanceyx  = sqrt( yx(:,1).^2 + yx(:,2).^2 );
            costhetayx  = sum( normal .* yx, 2 ) ./ distanceyx;
            
            usol2 = - log( distanceyx );
            dusol2 = - costhetayx ./ distanceyx;
            
            v2plot(it,ir) = sum( kernel .* ( 1- dusol2) .* mu .* Jacobian ) * dt  ...
                +  sum( kernel .* (mu - mu(it)) .* dusol2 .* Jacobian ) * dt  ...
                + mu(it) *  sum( kernel_dlp .* (usol2 - usol2(it)) .* Jacobian ) * dt ;
            
            %% METHOD 3: quadratic difference density subtraction
            coef3 = normal(it,1) .* rB(it,1) - normal(it,2) .* rB(it,2) ;
            usol3 = ( (rB(:,1).* rB(:,1)) - (rB(:,2).* rB(:,2))) / (2* coef3);
            dusol3 = (normal(:,1) .* rB(:,1) - normal(:,2) .* rB(:,2))/coef3;
            
            v3plot(it,ir) = sum( kernel .* (1 - dusol3) .* mu .* Jacobian ) * dt  ...
                +  sum( kernel .* (mu - mu(it)) .* dusol3 .* Jacobian ) * dt  ...
                + mu(it) *  sum( kernel_dlp .* (usol3 - usol3(it)) .* Jacobian ) * dt ;
            
            %% METHOD 4: quadratic product density subtraction
            coef4 = normal(it,1) .* (rB(it,2) - 5) + normal(it,2) .* (rB(it,1) -5) ;
            usol4 = (rB(:,1) -5) .* (rB(:,2)-5) / coef4;
            dusol4 = (normal(:,1) .* (rB(:,2) -5) + normal(:,2) .* (rB(:,1)-5))/coef4;
            
            v4plot(it,ir) = sum( kernel .* (1 - dusol4) .* mu .* Jacobian ) * dt  ...
                +  sum( kernel .* (mu - mu(it)) .* dusol4 .* Jacobian ) * dt  ...
                + mu(it) *  sum( kernel_dlp .* (usol4 - usol4(it)) .* Jacobian ) * dt ;
            
            %% METHOD 5: exponential density subtraction
            coef5 = exp(rB(it,1)) * (normal(it,1) .* sin(rB(it,2)) + normal(it,2) .* cos(rB(it,2))) ;
            usol5 = (  exp(rB(:,1)) .* sin(rB(:,2)) ) / (coef5);
            dusol5 = exp(rB(:,1)) .* (normal(:,1) .* sin(rB(:,2)) + normal(:,2) .* cos(rB(:,2))) /coef5;
            
            v5plot(it,ir) = sum( kernel .* (1 - dusol5) .* mu .* Jacobian ) * dt  ...
                +  sum( kernel .* (mu - mu(it)) .* dusol5 .* Jacobian ) * dt  ...
                + mu(it) *  sum( kernel_dlp .* (usol5 - usol5(it)) .* Jacobian ) * dt ;
            
            %% evaluate the exact solution
            
            wplot(it,ir) = ( xplot(it,ir) - x0 ) ...
                ./ ( ( xplot(it,ir) - x0 )^2 + ( yplot(it,ir) - y0 )^2 );
            
        end
        
    end
    %% report the maximum errors
    disp( ['   Maximum absolute error made by native PTR: ', ...
        num2str( max( max( abs( v0plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by linear density sub: ', ...
        num2str( max( max( abs( v1plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by Green function density sub: ', ...
        num2str( max( max( abs( v2plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    
    disp( ['   Maximum absolute error made by quadratic difference density sub: ', ...
        num2str( max( max( abs( v3plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by quadratic product density sub: ', ...
        num2str( max( max( abs( v4plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    disp( ['   Maximum absolute error made by exponential and sin density sub: ', ...
        num2str( max( max( abs( v5plot - wplot ) ) ) ) ] );
    
    disp( ' ' );
    
    
    %% impose periodicity
    
    mu(M+1)      = mu(1);
    xplot(M+1,:) = xplot(1,:);
    yplot(M+1,:) = yplot(1,:);
    v0plot(M+1,:) = v0plot(1,:);
    v1plot(M+1,:) = v1plot(1,:);
    v2plot(M+1,:) = v2plot(1,:);
    v3plot(M+1,:) = v3plot(1,:);
    v4plot(M+1,:) = v4plot(1,:);
    v5plot(M+1,:) = v5plot(1,:);
    wplot(M+1,:) = wplot(1,:);
    
    Err0(idx,:) =  ((abs( v0plot(it0,:) - wplot(it0,:))));
    Err1(idx,:) =  ((abs( v1plot(it0,:) - wplot(it0,:))));
    Err2(idx,:) =  ((abs( v2plot(it0,:) - wplot(it0,:))));
    Err3(idx,:) =  ((abs( v3plot(it0,:) - wplot(it0,:))));
    Err4(idx,:) =  ((abs( v4plot(it0,:) - wplot(it0,:))));
    Err5(idx,:) =  ((abs( v5plot(it0,:) - wplot(it0,:))));
end

%% Error slice

points = [M/8, M/2+1, M/2 + 30];
pts_names = ['A','B','C'];
for idx = 1: 3
    it0 = points(idx);
    f = figure;
    loglog( rgrid(1:end), abs( v0plot(it0,:) - wplot(it0,:) ), 'k--', ...
        rgrid(1:end), abs( v1plot(it0,:) - wplot(it0,:) ), '--x', ...
        rgrid(1:end), abs( v2plot(it0,:) - wplot(it0,:) ), '--o' ,...
        rgrid(1:end), abs( v3plot(it0,:) - wplot(it0,:) ), '--.',...
        rgrid(1:end), abs( v4plot(it0,:) - wplot(it0,:) ), '--+');
    xlabel( '$\ell$', 'Interpreter', 'LaTeX' );
    ylabel( 'Absolute error','Interpreter', 'LaTeX' );
    legend( 'V0', 'V1','V2', 'V3', 'V4','Interpreter', 'LaTeX','Location', 'northeastoutside');
    title(['Error at point ',num2str(pts_names(idx))],'Interpreter', 'LaTeX');
    xlim( [ 1.e-10 1] )
    ylim( [ 1e-15 1e0 ] )
    grid on;
    f.Position = [10 10 550 400];
end

%% Plot VS N
for ll = 1:7
    f = figure;
    loglog(MM, Err0(:,end-ll), 'k--', ...
        MM,Err1(:,end-ll), '--x',...
        MM, Err2(:,end-ll), '--o',...
        MM, Err3(:,end-ll), '--.',...
        MM, Err4(:,end-ll), '--+');
    xlabel( '$N$', 'Interpreter', 'LaTeX' );
    ylabel( 'Absolute error', 'Interpreter', 'LaTeX' );
    title(['Error at point A, $\ell$ = ', num2str(rgrid(end-ll))],'Interpreter', 'LaTeX');
    grid on;
    xlim([MM(1) MM(end)])
    ylim([1.e-10 1])
    legend( 'V0','V1','V2', 'V3', 'V4', 'Interpreter', 'LaTeX','Location', 'northeastoutside');
    f.Position = [10 10 550 400];
end

% %% Plot of the errors
% % METHOD 0
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v0plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% text(rB(12 ,1)-0.1,rB(12,2)-0.2,'A','FontSize',20,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
% text(rB(M/2+1,1)+0.1,rB(M/2+1,2),'B','FontSize',20 ,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% text(rB(M/2+30,1)+0.1,rB(M/2+30,2)+0.1,'C','FontSize',20 ,'Interpreter', 'LaTeX');
% title( 'Error V0' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% saveas(f,['SLP2D_V0_M',num2str(M),'.png']);
% 
% % METHOD 1
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v1plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% text(rB(12 ,1)-0.1,rB(12,2)-0.2,'A','FontSize',20,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
% text(rB(M/2+1,1)+0.1,rB(M/2+1,2),'B','FontSize',20 ,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% text(rB(M/2+30,1)+0.1,rB(M/2+30,2)+0.1,'C','FontSize',20 ,'Interpreter', 'LaTeX');
% title( 'Error V1' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% saveas(f,['SLP2D_V1_M',num2str(M),'.png']);
% % METHOD 2
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v2plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2,2),'xk','LineWidth',2);
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% title( 'Error V2' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% saveas(f,['SLP2D_V2_M',num2str(M),'.png']);
% 
% % METHOD 3
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v3plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% text(rB(12 ,1)-0.1,rB(12,2)-0.2,'A','FontSize',20,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
% text(rB(M/2+1,1)+0.1,rB(M/2+1,2),'B','FontSize',20 ,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% text(rB(M/2+30,1)+0.1,rB(M/2+30,2)+0.1,'C','FontSize',20 ,'Interpreter', 'LaTeX');
% title( 'Error V3' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% saveas(f,['SLP2D_V3_M',num2str(M),'.png']);
% 
% % METHOD 4
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v4plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% text(rB(12 ,1)-0.1,rB(12,2)-0.2,'A','FontSize',20,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
% text(rB(M/2+1,1)+0.1,rB(M/2+1,2),'B','FontSize',20 ,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% text(rB(M/2+30,1)+0.1,rB(M/2+30,2)+0.1,'C','FontSize',20 ,'Interpreter', 'LaTeX');
% title( 'Error V4' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% saveas(f,['SLP2D_V4_M',num2str(M),'.png']);
% 
% % METHOD 5
% f= figure;
% f.Position = [10 10 550 400];
% surf(xplot, yplot,log10( abs( v5plot - wplot ) ) );
% shading interp; view(2); colorbar; caxis([-16 1]);
% hold on;
% scatter(rB(12 ,1),rB(12,2),'xk','LineWidth',2);
% text(rB(12 ,1)-0.1,rB(12,2)-0.2,'A','FontSize',20,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +1  ,1),rB(M/2+1,2),'xk','LineWidth',2);
% text(rB(M/2+1,1)+0.1,rB(M/2+1,2),'B','FontSize',20 ,'Interpreter', 'LaTeX');
% hold on;
% scatter(rB(M / 2 +30  ,1),rB(M/2 +30,2),'xk','LineWidth',2);
% text(rB(M/2+30,1)+0.1,rB(M/2+30,2)+0.1,'C','FontSize',20 ,'Interpreter', 'LaTeX');
% title( 'Error V5' ,'Interpreter', 'LaTeX');
% axis tight;
% axis square;
% 
% 
