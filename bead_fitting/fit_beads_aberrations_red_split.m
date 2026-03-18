% Fit aberrations from beads and plot results

% Input: z-stack of beads

% - First the raw data is segmented in ROIs
% - Then, the aberrations of each ROIs (bead) is fitted.
% - The results are plotted.

clear all
close all

% Settings

% File that contains the raw data
path_input_data = 'C:\Users\Aleksandra\Libraries\vectorfit\data\beads_red_split\avg_pos3_red_split.tif';

% Folder where the results will be stored
path_output_data = 'C:\Users\Aleksandra\Libraries\vectorfit\data\beads_red_split\';
params = set_parameters_TetraSpeck_beads_100_red_split;

% Set up parallel pool
if isempty(gcp('nocreate'))
    number_of_parallel_workers = params.nr_of_parallel_workers;
    parpool('Threads', number_of_parallel_workers);
end

all_input_files = dir(path_input_data);
nr_of_files = size(all_input_files,1);
fprintf('Found %i files\n',nr_of_files);

for file_i=1:nr_of_files
    filename = all_input_files(file_i).name; 
    foldername = all_input_files(file_i).folder;
    path_input_data_full = fullfile(foldername, filename);    
    fprintf('Input bead data = %s\n',path_input_data_full);

    % Segmentation
    segmentation_beads_difference_uniform  

    % Aberration fitting   
    Ncfg_total = size(allspots,4);
    params.Ncfg_total = Ncfg_total;
    params.Ncfg_global = Ncfg_total;
    params.Ncfg = Ncfg_total;
    
    fprintf('\nPerform local fit with %i spots\n',Ncfg_total);
    
    max_local_iterations = params.max_local_iterations;
    numparams = params.numparams;
    
    theta_global = zeros(13,1);
    
    theta_local = initialvalues_phasor(allspots,roixy,theta_global,params);
    thetastore_local = zeros(numparams,Ncfg_total,max_local_iterations);
    meritstore = zeros(Ncfg_total,max_local_iterations);
    alambda_local = ones(Ncfg_total,1)*params.alambda0_local;
    alambdastore_local = zeros(Ncfg_total,params.max_local_iterations);
    Niters = zeros(Ncfg_total,1);
    flip_z_net = false;
    iiter_total = 1;
    
    [theta_local,thetastore_local,localizations,localizations_with_outliers,meritstore,alambda_local,mu,dmudtheta,Niters,outliers] = ...
    local_update(theta_local,thetastore_local,theta_global,meritstore,alambda_local,allspots,roixy,iiter_total,Niters,flip_z_net,framelist,ID,params);
    
    theta.local = theta_local;
    theta.global = theta_global;
    
    [~,filename,~] = fileparts(all_input_files(file_i).name); 
    path_output_data_full = strcat(path_output_data,'results_',filename,'.mat');
    save(path_output_data_full,'theta','thetastore_local','roixy','mu','allspots','Niters','outliers','params','-v7.3')
    fprintf('Fitting of %i beads finished\n',Ncfg_total);

end

%% Load and combine data from multiple files
path_load_data = strcat(path_output_data,'results_*.mat');

all_load_files = dir(path_load_data);
nr_of_files = size(all_load_files,1);
fprintf('Found %i files\n',nr_of_files);

theta_combined = [];
roixy_combined = [];
mu_combined = [];
allspots_combined = [];

for file_i=1:nr_of_files

    filename = all_load_files(file_i).name; 
    foldername = all_load_files(file_i).folder;
    path_load_data_full = fullfile(foldername, filename);    
    
    dataout = load(path_load_data_full,'theta','roixy','mu','allspots','outliers','params');
    
    theta_temp = dataout.theta.local;
    roixy_temp = dataout.roixy;
    mu_temp = dataout.mu;
    allspots_temp = dataout.allspots;
    outliers = dataout.outliers;
    no_outliers = setdiff(1:size(theta_temp,2),outliers);
    params = dataout.params;

    theta_combined = cat(2,theta_combined,theta_temp(:,no_outliers));
    roixy_combined = cat(2,roixy_combined,roixy_temp(:,no_outliers));
    mu_combined = cat(4,mu_combined,mu_temp(:,:,:,no_outliers));
    allspots_combined = cat(4,allspots_combined,allspots_temp(:,:,:,no_outliers));

end
fprintf('Finished loading bead data\n');

%% Remove beads with a large fiterror

errM = get_fiterror(mu_combined,allspots_combined,params);
idx_err = find(errM(1,:) > mean(errM(1,:)) + 2*std(errM(1,:)) | errM(3,:) > mean(errM(3,:)) + 2*std(errM(3,:)));
small_err = setdiff(1:size(roixy_combined,2),[idx_err]);

theta_combined = theta_combined(:,small_err);
roixy_combined = roixy_combined(:,small_err);
mu_combined = mu_combined(:,:,:,small_err);
allspots_combined = allspots_combined(:,:,:,small_err);


%% Visualize measured and modeled z-stack for one bead
ibead = 1; % choose bead number
dipshow(cat(2,allspots(:,:,:,ibead),mu(:,:,:,ibead)),'lin')
diptruesize(1500)
colormap parula

%% Plot sample tilt (z-value of beads)
pixelsizex = params.pixelsizex;
pixelsizey = params.pixelsizey;
xsize = pixelsizex*params.imgSizeX;
ysize = pixelsizey*params.imgSizeY;

[Xq,Yq] = meshgrid(1:xsize*1e-3,1:ysize*1e-3);

% Plot sample tilt
Zq = griddata(roixy_combined(1,:)*pixelsizex*1e-3,roixy_combined(2,:)*pixelsizey*1e-3,theta_combined(3,:),Xq,Yq,"cubic");
figure
surf(Xq,Yq,Zq,'FaceAlpha',0.7,'EdgeColor','none');
hold on
plot3(roixy_combined(1,:)*pixelsizex*1e-3,roixy_combined(2,:)*pixelsizey*1e-3,theta_combined(3,:),"o",'LineWidth',2);
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (nm)')
title('Sample tilt')
ax = gca;
ax.FontSize = 32;

%% Plot aberrations mean and std
zernike_coeffs = theta_combined(6:end,:)*1e3/params.lambda;%% Plot aberration surfaces
figure
hold on; box on;
numzers = params.numparams-5;
plot(0:numzers+1,zeros(1,numzers+2),'-','Color',[.85 .85 .85],'LineWidth',0.5)
orders = params.aberrations(:,1:2);
allxticks = 1:numzers;
allxticklabels = cell(numzers,1);
for jzer = 1:numzers
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end
plot(1:numzers,zernike_coeffs,'k-*','MarkerSize',5)
xticks(allxticks)
xtickangle(25)
xticklabels(allxticklabels)
xlim([0 numzers+1])
xlabel('zernike mode (n,m)');
ylabel('m\lambda');
title('Aberrations per bead')

% Plot mean and standard deviation
mean_aber = mean(zernike_coeffs,2);
std_aber = std(zernike_coeffs,0,2);
figure
hold on; box on;
numzers = params.numparams-5;
plot(0:numzers+1,zeros(1,numzers+2),'-','Color',[.85 .85 .85],'LineWidth',0.5)
errorbar(1:numzers,mean_aber,std_aber,'Color','black','LineWidth',1)
plot(1:numzers,mean_aber,'k-*','MarkerSize',10)
xticks(allxticks)
xtickangle(25)
xticklabels(allxticklabels)
xlim([0 numzers+1])
xlabel('zernike mode (n,m)');
ylabel('m\lambda');
title('Mean and standard deviation aberrations')

%% Plot aberrations across the FOV
numzers = params.numparams-5;
orders = params.aberrations(:,1:2);
allxticks = 1:numzers;
allxticklabels = cell(numzers,1);
for jzer = 1:numzers
    allxticklabels{jzer} = strcat(num2str(orders(jzer,1)),',',num2str(orders(jzer,2)));
end

pixelsizex = params.pixelsizex;
pixelsizey = params.pixelsizey;
xsize = pixelsizex*params.imgSizeX;
ysize = pixelsizey*params.imgSizeY;

[Xq,Yq] = meshgrid(1:xsize*1e-3,1:ysize*1e-3);

zernike_coeffs = theta_combined(6:end,:)*1e3/params.lambda;%% Plot aberration surfaces

[~,FOV_coordinates] = get_fov_coordinates(roixy_combined,theta_combined(1,:),theta_combined(2,:),params);
X = FOV_coordinates(1,:);
Y = FOV_coordinates(2,:);

zmin = min(zernike_coeffs,[],2);
zmax = max(zernike_coeffs,[],2);
zmin = zmin - 0.5*(zmax - zmin);
zmax = zmax + 0.5*(zmax - zmin);

figure('units','normalized','outerposition',[0 0 1 1])
for izer = 1:size(params.aberrations,1)
    hsub = subplot(2,3,izer);
    V = zernike_coeffs(izer,:);
    Vq = griddata(X,Y,V,Xq,Yq,"cubic");

    surf(Xq,Yq,Vq,'FaceAlpha',0.7)
    shading interp
    hold on
    plot3(X,Y,V,"o",'LineWidth',3,'MarkerSize',10)
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('m\lambda')
    zlim([zmin(izer) zmax(izer)])
    title(strcat('A(',allxticklabels{izer},')'))
    ax = gca;
    ax.FontSize = 16;

end
