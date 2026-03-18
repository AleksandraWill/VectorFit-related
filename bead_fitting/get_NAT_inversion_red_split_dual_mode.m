% Interpolate aberration surfaces via inversion of NAT

clear all
close all

%%
path_output_data = 'C:\Users\Aleksandra\Libraries\vectorfit\data\beads_red_blue_split_red_dual\results_*';

all_output_files = dir(path_output_data);
nr_of_outfiles = size(all_output_files,1);
fprintf('Found %i files\n',nr_of_outfiles);

theta_combined = [];
roixy_combined = [];
mu_combined = [];
allspots_combined = [];

for file_i=1:nr_of_outfiles
%for file_i=6:9

    filename = all_output_files(file_i).name
    foldername = all_output_files(file_i).folder;
    path_output_data_full = fullfile(foldername, filename);    
    
    dataout = load(path_output_data_full,'theta','roixy','mu','allspots','outliers','params');
    
    theta_temp = dataout.theta.local;
    roixy_temp = dataout.roixy;
    mu_temp = dataout.mu;
    allspots_temp = dataout.allspots;
    outliers = dataout.outliers;
    %outliers = [];
    no_outliers = setdiff(1:size(theta_temp,2),outliers);
    params = dataout.params;

    theta_combined = cat(2,theta_combined,theta_temp(:,no_outliers));
    roixy_combined = cat(2,roixy_combined,roixy_temp(:,no_outliers));
    mu_combined = cat(4,mu_combined,mu_temp(:,:,:,no_outliers));
    allspots_combined = cat(4,allspots_combined,allspots_temp(:,:,:,no_outliers));

end

% compatibility with old pixelsize convention
if exist('params', 'var') && isfield(params, 'pixelsize')
    params.pixelsizex = params.pixelsize;
    params.pixelsizey = params.pixelsize;
end

errM = get_fiterror(mu_combined,allspots_combined,params);

idx_err = find(errM(1,:)>50 | errM(2,:)>50 | errM(3,:)>5e4); % Lidke 2D

small_err = setdiff(1:size(roixy_combined,2),idx_err);

theta_combined = theta_combined(:,small_err);
roixy_combined = roixy_combined(:,small_err);
mu_combined = mu_combined(:,:,:,small_err);
allspots_combined = allspots_combined(:,:,:,small_err);
%%
params.aberrations = [...
     2     0     0
     2    -2     0
     2     2     0
     3    -1     0
     3     1     0
     4     0     0];

params.NATgammas = {{[0 0 2 0 1.0],0.0},...
                     {[1 0 2 0 1.0],0.0},...
                     {[0 1 2 0 1.0],0.0},...
                     {[2 0 2 0, 1.0],[0 2 2 0, 1.0],0.0},...
                     {[0 0 2 2 1.0],0.0},...
                     {[0 0 2 -2 1.0],0.0},...
                     {[1 0 2 2 1.0],[0 1 2 -2 1.0],0.0},...
                     {[0 1 2 2 -1.0],[1 0 2 -2 1.0],0.0},...
                     {[2 0 2 2 1/sqrt(5)],[0 2 2 2 -1/sqrt(5)],[1 1 2 -2 1.0],0.0},...
                     {[0 0 3 1 1.0],0.0},...
                     {[0 0 3 -1 1.0],0.0},...
                     {[1 0 3 1 1.0],[0 1 3 -1 1.0],0.0},...
                     {[0 0 4 0 1.0],0.0}};

params.numgammas = length(params.NATgammas);
params.multiplex_factors = zeros(params.numgammas,1); 
params.numNATparvecs = 0;
for jgam = 1:params.numgammas
  RR = params.NATgammas{jgam};
  multiples = length(RR)-1;
  for jm = 1:multiples
    NATvec = RR{jm};
    params.multiplex_factors(jgam) = params.multiplex_factors(jgam) + NATvec(5)^2;
  end
  params.numNATparvecs = params.numNATparvecs + multiples;
end

params.orders2D = zeros(params.numNATparvecs,2);
    jv = 1; % counter to step through list of all different NATvec terms
    for jgam = 1:params.numgammas
      RR = params.NATgammas{jgam};
      multiples = length(RR)-1;
      for jm = 1:multiples
        NATvec = RR{jm};
        params.orders2D(jv,1) = NATvec(1);
        params.orders2D(jv,2) = NATvec(2);
        jv = jv+1;
      end
    end

params.legendre_normfac = sqrt((1+2*params.orders2D(:,1)).*(1+2*params.orders2D(:,2)));


%%
Ncfg = size(roixy_combined,2);
allzernikes_measured = [zeros(1,Ncfg); theta_combined(6:end,:)];
fov_coordinates = get_fov_coordinates(roixy_combined,theta_combined(1,:),theta_combined(2,:),params);
xn = fov_coordinates(1,:);
yn = fov_coordinates(2,:);

%%
xgrid = linspace(-1,1,50);
ygrid = linspace(-1,1,50);
[Xgrid,Ygrid] = meshgrid(xgrid,ygrid);

gammas = invertNAT(xn,yn,allzernikes_measured,params);
allzernikes_fitted = 1e3*get_zernike_coefficients(xn,yn,gammas,params)'/params.lambda;
zernike_surface_fitted = 1e3*get_zernike_coefficients(Xgrid,Ygrid,gammas,params)/params.lambda;

%%

labels = ["A(2,0)" "A(2,-2)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];
figure
for izer = 2:size(params.aberrations,1)
    hsub = subplot(2,3,izer);
    V = 1e3*allzernikes_measured(izer,:)/params.lambda;
    Vgrid = griddata(xn,yn,V,Xgrid,Ygrid,"cubic");

    surf(Xgrid,Ygrid,Vgrid,'FaceAlpha',0.7)
    shading interp
    hold on
    plot3(xn,yn,V,"o",'LineWidth',2,'MarkerSize',5)
    hold on
    surf(Xgrid,Ygrid,zernike_surface_fitted(:,:,izer),'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#F80')
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('m\lambda')
    %zlim([-40 60])
    title(labels(izer))
    ax = gca;
    ax.FontSize = 16;

end

legend('Bead aberrations cubic interpolation','Beads','NAT surfaces from beads','Location','east','Fontsize',16)


