 % Matlab segmentation for single emitters in SMLM data 
% Based on Huanget al. (2011)
% 
% Fang Huang, Samantha L. Schwartz, Jason M. Byars, and Keith A. Lidke, 
% "Simultaneous multiple-emitter fitting for single molecule super-resolution 
% imaging," Biomed. Opt. Express 2, 1377-1393 (2011)
% -------------------------------------------------------------------------
% set parameters

% uncomment to run this file separately
%clear all
%close all
%params = set_parameters_zstack_bead_Cali1031_100bead_Oil100x;
%path_input_data = 'N:\tnw\IST\QI\users\idroste\Data\beads_example_for_yutong\Nanoruler\Raw data\CalibBead\Cali1031_100bead_Oil100x_Col0-17_ex488_P50_400ms_EM100_Col170-Cntr.tif';
%path_output_data = 'N:\tnw\IST\QI\users\idroste\Data\beads_example_for_yutong\Nanoruler\Raw data\CalibBead\Output\Cali1031_100bead_Oil100x_Col0-17_ex488_P50_400ms_EM100_Col170-Cntr_IsabelAnalysis\segmentation\segmentation_file_';

% Use when a gain and offset image are provided.
if params.gain_image
    %path_gain_image = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Hertenlab\raw_data\CameraCalibrationMap\Corrected_Gain2_Balanced.tif';
    path_gain_image = '/home/idroste/Desktop/TUDelft/Data/Hertenlab/raw_data/CameraCalibrationMap/Corrected_Gain2_Balanced.tif';
    fprintf('Use provided gain image %s\n',path_gain_image)
    gain = imread(path_gain_image);
else
    gain = params.gain;
end
if params.offset_image
    %path_offset_image = 'U:\SMLMaberrations\NATaberrations\ExperimentalData\Hertenlab\raw_data\CameraCalibrationMap\Offset_Balanced_ADU.tif';
    path_offset_image = '/home/idroste/Desktop/TUDelft/Data/Hertenlab/raw_data/CameraCalibrationMap/Offset_Balanced_ADU.tif';
    fprintf('Use provided offset image %s\n',path_offset_image)
    offset = imread(path_offset_image);
else
    offset = params.offset;
end

sigma_psf = params.lambda/(4*params.NA*sqrt(params.pixelsizex*params.pixelsizey));
segmentation_ROI_thr = params.segmentation_ROI_thr;
unif1=round_odd(2*sigma_psf+1);
unif2=round_odd(2*unif1);
maxf3=round(5*sigma_psf);

roisize = params.Mx; % has to be odd
ROI_pad = (roisize-1)/2;
imgsizeX = params.imgSizeX;
imgsizeY = params.imgSizeY;
Mz = params.Mz;
%maxspots = 250E3;
%maxframe = 5E3;

fprintf('\nStart segmentation')

tic

fprintf('\nProcess file: %s \n', path_input_data_full)

allroi_centers = zeros(imgsizeY,imgsizeX);

%switch params.datatype
%    case('tif')
%    InfoImage =  imfinfo(path_input_data_full);
%    Nframes = length(InfoImage);
%    
%    fprintf('Loading %i frames\n',Nframes);
%    allframes = double(LoadTiff(path_input_data_full, Nframes));
%    fprintf('Finished loading %i frames\n',Nframes)
%    
%    case('mat')
%        fprintf('Loading data\n');
%        input_data = load(path_input_data_full);
%        allframes = input_data.ZStack;
%        Nframes = size(allframes,3);
%        fprintf('Finished loading %i frames\n',Nframes)
%    otherwise
%        error('params.datatype not correct')
%end

allframes = (allframes-offset)./gain;
allframes(allframes<=0) = 1e-3;
Nframes = size(allframes,3);
frame = allframes(:,:,round(size(allframes,3)/2));

box_width_x = size(frame,1);
box_width_y = size(frame,2);
A1 = imboxfilt(frame,unif1) - imboxfilt(frame,unif2);

A2 = imdilate(A1,ones(maxf3,maxf3));
A3 = (A1 == A2) & (A1 > segmentation_ROI_thr);
ind_uf = find(A3);

% Find roixy
nrois_all = length(ind_uf);
roixy_unfiltered = zeros(2,nrois_all);
for ii = 1:nrois_all
    [idx_uf, idy_uf] = ind2sub(size(frame), ind_uf(ii));
    roixy_unfiltered(:,ii) = [idx_uf idy_uf]';
end

% Remove overlapping ROIs
mindist = ceil(roisize/2);
distances = squareform(pdist(roixy_unfiltered'));
[distances_small_x,distances_small_y] = ind2sub([nrois_all nrois_all],find(distances<=mindist & distances>0));

remove_indices = [-1];
for i = 1:numel(distances_small_y)
    idy = distances_small_y(i);
    idx = distances_small_x(i);
    if ~ismember(remove_indices,idy)
        remove_indices = cat(1,remove_indices,idx);
    end
end
remove_indices = unique(remove_indices(2:end));
indices_keep = setdiff(1:nrois_all,remove_indices);
ind = ind_uf(indices_keep);
roixy = roixy_unfiltered(:,indices_keep);

% Cut out filtered ROIs
nrois = length(ind);
allspots = zeros(roisize,roisize,Mz,nrois);
framelist = ones(1,nrois);
for ii = 1:nrois
    [idx, idy] = ind2sub(size(frame), ind(ii));
    roixy(:,ii) = [idx idy]';
    [idx_ROI, idy_ROI] = ROI_coords(idx, idy, box_width_x, box_width_y, ROI_pad);

    ROI = allframes(min(idx_ROI) : max(idx_ROI), min(idy_ROI) : max(idy_ROI),:);
    allspots(:,:,:,ii) = ROI;
end

Ncfg_total = nrois;
ID = 1:Ncfg_total;
params.Ncfg_total = Ncfg_total;
params.Mz = Nframes;
fprintf('Found %i spots\n',Ncfg_total);
%path_output_data_full = strcat(path_output_data,num2str(file_i),'_segmentation_.mat');
%save(path_output_data_full,'allspots','roixy','framelist','ID','params');
toc

fprintf('\nFinished segementation\n');

%% Plots segmentation

% Plot all ROI centers
allroi_centers = zeros(imgsizeY,imgsizeX);
for i = 1:size(roixy,2)
    allroi_centers(roixy(1,i),roixy(2,i)) = 1;
end

dipshow(allroi_centers,'lin')
diptruesize(5e4/params.imgSizeX)

dipshow(frame,'lin')
diptruesize(5e4/params.imgSizeX)

figure('Name','Use this figure to check the threshold')
surf(A1);
shading interp;

dipshow(allspots(:,:,round(size(allframes,3)/2),:),'lin')
diptruesize(4e4/params.Mx)
colormap parula

% Xim = 1e-3*params.pixelsizex*[1:params.imgSizeX];
% Yim = 1e-3*params.pixelsizey*[1:params.imgSizeY];
% frame_img = mat2im(frame);
% frame_gauss = gaussf(frame_img,15);
% frame_gauss_mat = im2mat(frame_gauss);
% figure('Name','Frame illumination')
% surf(Xim,Yim,frame_gauss_mat,'FaceAlpha',1.0);
% shading interp;
% xlabel('x (\mum)');
% ylabel('y (\mum)');
% zlabel('Intensity');
% ax=gca;
% ax.FontSize = 18;

function S = round_odd(S)
    % round to nearest odd integer.
    idx = mod(S,2)<1;
    S = floor(S);
    S(idx) = S(idx)+1;
end
