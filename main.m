% Sample code from "From sparse data to high-resolution fields: ensemble 
% particle modes as a basis for high-resolution flow characterization".
%
% J. Cortina Fernández, C. Sanmiguel Vila, A. Ianiro, S. Discetti.

clc, clear, close all

%% Parameters
% Tunable parameters
Np    = 800;    % Number of particles per bin
r     = 8;      % Number of modes of the high-resolution reconstruction
r_PIV = 8;      % Number of modes of the PIV reconstruction

% Parameters that depend on the INPUT data
Nppp = 0.02;    % Particles per pixel
dx   = 4;       % Grid size [px]
Nm   = 150;     % Number of high-resolution POD modes to be computed
Nt   = 1510;    % Number of images used for reconstruction
sten = 2;       % Stencil size used in the spatial filter [px]

% Bin size
b = sqrt(Np/(Nt*Nppp)); % [px]

% Resolution, number of pixels per diameter
res = 25;   % [px/D]

% Image resolution
nx = 225;   % [px]
ny = 99;    % [px]

% Cylinder location in the image [x,y]
cyl_pos = [25, 50]; % [px]

%% Load data
% INPUT data
PIV = load('.\INPUT\PIV.mat');  % PIV data
PTV = load('.\INPUT\PTV.mat');  % PTV data

% Auxiliary functions
addpath('.\functions\')
addpath('.\functions\kdtree\')

%% DEPTV algorithm
% Temporal basis: Perform the Proper Orthogonal Decomposition (POD) of the
% fluctuating component of the PIV data to obtain the temporal basis of the
% DEPTV algorithm
    % Time-averaged velocity field
    PIV.Um = mean(PIV.U,3);
    PIV.Vm = mean(PIV.V,3);
    
    % Fluctuating velocity field
    PIV.u = PIV.U - PIV.Um;
    PIV.v = PIV.V - PIV.Vm;
    
    % Reshape matrix from 3D to 2D where: rows = time & cols = gridpoints
    PIV.u = reshape(PIV.u,[numel(PIV.X),Nt])';
    PIV.v = reshape(PIV.v,[numel(PIV.X),Nt])';
    
    % Perform POD on snapshot matrix of the fluctuating PIV data
    [PIV.psi,PIV.sigma,PIV.phi] = svd([PIV.u,PIV.v],'econ');

% Spatial basis: Calculate the ensemble-particle modes
    % Define DEPTV grid
    x     = dx:dx:(nx-dx);                                          % [px]
    y     = dx:dx:(ny-dx);                                          % [px]
    [X,Y] = meshgrid((x - cyl_pos(1))/res,(y - cyl_pos(2))/res);    % [D]
    
    % Calculate time-averaged velocity field
    fprintf('Calculating time-averaged PTV field:\n')
    
        % Group all particles into "_all" variables
        PTV.X_all = PTV.X{1};
        PTV.Y_all = PTV.Y{1};
        PTV.U_all = PTV.U{1};
        PTV.V_all = PTV.V{1};
        for i = 2:Nt
            if mod(i,floor(Nt/10)) == 0
                fprintf('    %d%%\n',i/Nt*100)
            end
            PTV.X_all = [PTV.X_all;PTV.X{i}];
            PTV.Y_all = [PTV.Y_all;PTV.Y{i}];
            PTV.U_all = [PTV.U_all;PTV.U{i}];
            PTV.V_all = [PTV.V_all;PTV.V{i}];
        end

        % Build k-d tree
        tree = kdtree_build([PTV.X_all(:),PTV.Y_all(:)]); 

        % Find all particles inside bin for each grid point
        PTV.Um = zeros(1,numel(X));
        PTV.Vm = PTV.Um;
        for i = 1:numel(X)  % For each grid point
            % Find particles inside bin
            range = [X(i) + b*[-0.5,0.5]/res; Y(i) + b*[-0.5,0.5]/res]; % [D]
            ij    = kdtree_range_query(tree,range);

            % Calculate velocity mean of particles within each bin
            PTV.Um(i) = mean(PTV.U_all(ij));
            PTV.Vm(i) = mean(PTV.V_all(ij));
        end

        % Delete tree
        kdtree_delete(tree);
    
    % Fluctuating velocity: substract time-averaged field interpolated to particles' locations
    PTV.u = cell(size(PTV.U));
    PTV.v = cell(size(PTV.V));
    for i = 1:Nt
        PTV.u{i} = PTV.U{i} - interp2(X,Y,reshape(PTV.Um,size(X)),PTV.X{i},PTV.Y{i},'spline');
        PTV.v{i} = PTV.V{i} - interp2(X,Y,reshape(PTV.Vm,size(X)),PTV.X{i},PTV.Y{i},'spline');
    end
    
    % Velocity projection into PIV temporal basis
    fprintf('Calculating projections:\n')
    sigphiu = zeros(Nm,numel(X));   % sigma*phi matrix of streamwise coordinate
    sigphiv = sigphiu;              % sigma*phi matrix of normal coordinate
    Nocc    = zeros(1,numel(X));    % Number of occurrences per bin
    for k = 1:Nt    % For each snapshot
        if mod(k,floor(Nt/10)) == 0
            fprintf('    %d%%\n',k/Nt*100)
        end
        % Build k-d tree of particles at snapshot 'k'
        tree = kdtree_build([PTV.X{k},PTV.Y{k}]); 

        for i = 1:numel(X)  % For each bin
            % Find particles inside bin of grid point 'i'
            range = [X(i) + b*[-0.5,0.5]/res; Y(i) + b*[-0.5,0.5]/res];
            ij    = kdtree_range_query(tree,range);
            
            % For each bin, if there are particles
            if ~isempty(ij)
                % Build the high-resolution ensemble-particle modes
                sigphiu(1:Nm,i) = sigphiu(1:Nm,i) + repmat(PIV.psi(k,1:Nm)',[1,numel(ij)])*PTV.u{k}(ij);
                sigphiv(1:Nm,i) = sigphiv(1:Nm,i) + repmat(PIV.psi(k,1:Nm)',[1,numel(ij)])*PTV.v{k}(ij);
                
                % Update occurrences
                Nocc(i) = Nocc(i) + numel(ij);    
            end
        end
        
        % Delete tree
        kdtree_delete(tree);
    end
    
    % Scale spatial modes using the number of occurrences per bin
    sigphiu = sigphiu./repmat(Nocc,[Nm,1]).*Nt;
    sigphiv = sigphiv./repmat(Nocc,[Nm,1]).*Nt;
            
% Apply Savitzky-Golay filter in ensemble-particle modes
    % Find grid points inside cylinder
    theta = 2*pi*(0:100)/100;               % Define circumference
    xCyl  = 0.5*cos(theta);                 % x coordinate
    yCyl  = 0.5*sin(theta);                 % y coordinate
    in    = inpolygon(X(:),Y(:),xCyl,yCyl); % Check for grid points inside circles 

    % Memory pre-allocation of filtered modes
    sigphiuF = sigphiu;
    sigphivF = sigphiv;

    % Filter each mode
    fprintf('Filtering modes\n')
    for i=1:Nm
        Vf = polyfilt(reshape(sigphiu(i,:),size(X)),sten,reshape(in,size(X)));
        sigphiuF(i,:) = Vf(:);
        Vf = polyfilt(reshape(sigphiv(i,:),size(X)),sten,reshape(in,size(X)));
        sigphivF(i,:) = Vf(:);
    end

% Instantaneous velocity field reconstruction with 'r_PIV' or 'r' number of modes
    % PIV fluctuating velocity field
    u_PIV = PIV.psi(:,1:r_PIV)*PIV.sigma(1:r_PIV,1:r_PIV)*PIV.phi(1:(size(PIV.phi,1)/2),1:r_PIV)';
    v_PIV = PIV.psi(:,1:r_PIV)*PIV.sigma(1:r_PIV,1:r_PIV)*PIV.phi((size(PIV.phi,1)/2+1):end,1:r_PIV)';
    
    % DEPTV fluctuating velocity field
    u_DEPTV = PIV.psi(:,1:r)*sigphiu(1:r,:);
    v_DEPTV = PIV.psi(:,1:r)*sigphiv(1:r,:);
    
    % Filtered DEPTV fluctuating velocity field
    u_DEPTVf = PIV.psi(:,1:r)*sigphiuF(1:r,:);
    v_DEPTVf = PIV.psi(:,1:r)*sigphivF(1:r,:);

%% Figures
fprintf('Creating figures\n')

% Singular value distribution
    % Calculate norm of the unfiltered and filtered ensemble-particle modes
    sigma_DEPTV = zeros(Nm,1);
    sigma_DEPTVf = zeros(Nm,1);
    for i = 1:Nm    % For each mode
        sigma_DEPTV(i)  = norm([sigphiu(i,:),sigphiv(i,:)]);
        sigma_DEPTVf(i) = norm([sigphiuF(i,:),sigphivF(i,:)]);
    end
    
    % Plot
    lw = 1.5;   % Linewidth
    figure
    loglog((diag(PIV.sigma).^2)./sum(sigma_DEPTV.^2)*(numel(X)/numel(PIV.X)),'o-b','LineWidth',lw)
    hold on
    loglog((sigma_DEPTV).^2/sum(sigma_DEPTV.^2),'^-','color',[0 0.5 0],'LineWidth',lw)
    loglog((sigma_DEPTVf).^2/sum(sigma_DEPTV.^2),'s-r','LineWidth',lw)
    grid on
    xlim([1,150])
    xticks([1 5 10 20 50 100 150])
    ylim([1e-5 1])
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1])
    legend('PIV','Raw DEPTV','Filtered DEPTV','interpreter','latex','location','southeast')
    set(gca,'FontSize',12,'ticklabelinterpreter','latex');
    xlabel('Mode number, $i$','FontSize',12,'interpreter','latex')
    ylabel('$\hat{\sigma}_i^2$','FontSize',14,'interpreter','latex')
    axes('Position',[.52 .55 .35 .33])
    box on
    semilogx((diag(PIV.sigma).^2)./sum(sigma_DEPTV.^2)*(numel(X)/numel(PIV.X)),'o-b','LineWidth',lw)
    hold on
    semilogx((sigma_DEPTV).^2/sum(sigma_DEPTV.^2),'^-','color',[0 0.5 0],'LineWidth',lw)
    semilogx((sigma_DEPTVf).^2/sum(sigma_DEPTV.^2),'s-r','LineWidth',lw)
    grid on
    xlim([1,8])
    xticks([1 2 4 6 8])
    yticks(0:0.1:0.5)
    ylim([0,0.5])
    set(gca,'FontSize',12,'ticklabelinterpreter','latex');
    set(gcf,'units','points','Position',[100 100 380 300]);

% Spatial modes 2, 4 and 6 projected on the PIV temporal basis
    % Modes being represented (Choose 3 modes)
    modes  = [2,4,6];
    
	% Spacing settings
    gap    = [0.07 0.03]; % [vertical,horizontal] defining the gap between neighbouring axes
    marg_h = [0.07 0.04]; % [lower uppper] margins in height in normalized units (0...1)
    marg_w = [0.05 0.07]; % [left right] margins in width in normalized units (0...1)
    
    % Plot
    figure
    for i = 1:length(modes)
        % PIV
        subtightplot(3,3,i,gap,marg_h,marg_w)
        imagesc(PIV.X(1,:),PIV.Y(:,1),reshape(PIV.phi(1:(size(PIV.phi,1)/2),modes(i)),size(PIV.X))*sqrt(numel(PIV.X)))
        title(['PIV mode ',num2str(modes(i))],'FontSize',12,'interpreter','latex')
        
        % Unfiltered DEPTV
        subtightplot(3,3,3+i,gap,marg_h,marg_w)
        imagesc(X(1,:),Y(:,1),reshape(sigphiu(modes(i),:)/sigma_DEPTVf(modes(i)),size(X))*sqrt(numel(X)))
        title(['Raw DEPTV mode ',num2str(modes(i))],'FontSize',12,'interpreter','latex')
        
        % Filtered DEPTV
        subtightplot(3,3,6+i,gap,marg_h,marg_w)
        imagesc(X(1,:),Y(:,1),reshape(sigphiuF(modes(i),:)/sigma_DEPTVf(modes(i)),size(X))*sqrt(numel(X)))
        title(['Filtered DEPTV mode ',num2str(modes(i))],'FontSize',12,'interpreter','latex')
    end
    
    for i = 1:9
        subtightplot(3,3,i,gap,marg_h,marg_w)
        hold on
        colormap jet
        caxis(1.5*[-1,1])
        daspect([1,1,1])
        fill(xCyl,yCyl,[1 1 1])             % place cylinder
        plot(xCyl,yCyl,'k','LineWidth',1.2) % cylinder boundary
        axis([-0.85 7.85 -1.8 1.8])
        ax = gca;
        ax.YDir = 'normal';
        set(gca,'FontSize',11,'ticklabelinterpreter','latex');
        if i == 7 || i == 8 || i == 9
            xlabel('$x/D$','FontSize',12,'interpreter','latex')
        end
        if i == 1 || i == 4 || i == 7
            ylabel('$y/D$','FontSize',12,'interpreter','latex')
        end
    end
    hcb = colorbar('location','Manual','position',[0.95 0.09 0.015 0.855],'FontSize',11,'ticklabelinterpreter','latex');
    set(gcf,'renderer','painters','units','points','Position',[100 100 700 370]);
    title(hcb,'$\phi_u \sqrt{2 p}$','Interpreter','latex','FontSize',14)

% Instant fluid field
    % Time instant being represented (Choose between 1 and Nt)
    t = 1;

    % Spacing settings
    gap    = [0.07 0.06];   % [vertical,horizontal] defining the gap between neighbouring axes
    marg_h = [0.07 0.04];   % [lower uppper] margins in height in normalized units (0...1)
    marg_w = [0.10 0.15];   % [left right] margins in width in normalized units (0...1)

    % Plot
    figure
        % PIV
        subtightplot(3,1,1,gap,marg_h,marg_w)
        imagesc(PIV.X(1,:),PIV.Y(:,1),reshape(u_PIV(t,:),size(PIV.X)))
        hold on
        quiver(PIV.X,PIV.Y,reshape(u_PIV(t,:),size(PIV.X)),reshape(v_PIV(t,:),size(PIV.X)),'k','Autoscale',2)
        title('PIV','FontSize',11,'interpreter','latex')

        % Raw DEPTV
        subtightplot(3,1,2,gap,marg_h,marg_w)
        imagesc(X(1,:),Y(:,1),reshape(u_DEPTV(t,:),size(X)))
        hold on
        quiver(X,Y,reshape(u_DEPTV(t,:),size(X)),reshape(v_DEPTV(t,:),size(X)),'k','Autoscale',2)
        title('Raw DEPTV','FontSize',11,'interpreter','latex')

        % Filtered DEPTV
        subtightplot(3,1,3,gap,marg_h,marg_w)
        imagesc(X(1,:),Y(:,1),reshape(u_DEPTVf(t,:),size(X)))
        hold on
        quiver(X,Y,reshape(u_DEPTVf(t,:),size(X)),reshape(v_DEPTVf(t,:),size(X)),'k','Autoscale',2)
        title('Filtered DEPTV','FontSize',11,'interpreter','latex')
        xlabel('$x/D$','FontSize',11,'interpreter','latex')

    for i = 1:3
        subtightplot(3,1,i,gap,marg_h,marg_w)
        hold on
        colormap jet
        caxis(0.5*[-1,1])
        daspect([1,1,1])
        fill(xCyl,yCyl,[1 1 1])             % place cylinder
        plot(xCyl,yCyl,'k','LineWidth',1.2) % cylinder boundary
        ax = gca;
        axis([-0.85 7.85 -1.8 1.8])
        ax.YDir = 'normal';
        set(gca,'FontSize',11,'ticklabelinterpreter','latex');
        ylabel('$y/D$','FontSize',11,'interpreter','latex')
    end
    hcb=colorbar('location','Manual','position',[0.88 0.079 0.03 0.875],'FontSize',11,'ticklabelinterpreter','latex');
    set(gcf,'renderer','painters','units','points','Position',[100 100 300 400]);
    title(hcb,'$u/U_\infty$','Interpreter','latex','FontSize',14)

