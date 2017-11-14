%HW5
%GB comments
1a 100
1b 100
1c 100
1d 100
2yeast: 70 output image is not a nuclear/cell mask but an edge mask. You need to tell your script that the cellular space is not the same as the background. 
2worm: 100 invert the image
2bacteria: 95 This is pretty good. Instead of using imdilate, you could have used imfill. It would have filled your holes while also not expanding the edges of the mask. When you dilated, several masks end up touching. Alternatively, you could have followed up with imerode to fix this issue. 
2phase: 100 
Overall: 96


% Note. You can use the code readIlastikFile.m provided in the repository to read the output from
% ilastik into MATLAB.

%% Problem 1. Starting with Ilastik

% Part 1. Use Ilastik to perform a segmentation of the image stemcells.tif
% in this folder. Be conservative about what you call background - i.e.
% don't mark something as background unless you are sure it is background.
% Output your mask into your repository. What is the main problem with your segmentation?  

% part 1 and part 2 will be answered together.

% Part 2. Read you segmentation mask from Part 1 into MATLAB and use
% whatever methods you can to try to improve it. 

% segmentation data from the Ilastik analysis from part 1
data_seg_p1 = h5read('Segmentation_Nuclei_part1.h5', '/exported_data');
data_seg_p1 = squeeze(data_seg_p1);
imshow(data_seg_p1, []);

% uncertainty data from the Ilastik analysis from part 1
data_unc_p1 = h5read('uncertainty_part1.h5', '/exported_data');
data_unc_p1 = squeeze(data_unc_p1);
imshow(data_unc_p1, []);

imshowpair(data_seg_p1, data_unc_p1);

% visualizing fused candidates
BW = im2bw(data_seg_p1, graythresh(data_seg_p1));
CC = bwconncomp(BW);
stats = regionprops(CC, 'Area');
area = [stats.Area];
fusedCandidates = area > mean(area) + std(area);
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:});
fusedMask = false(size(data_seg_p1));
fusedMask(sublist) = 1;
imshow(fusedMask, 'InitialMagnification', 'fit');

t_p1 = cat(3, fusedMask, BW, zeros(size(BW)));
imshow(t_p1);

% As we can see from the data, using a very conservative analysis for
% background in Ilastik a lot of uncertainty is generated in the file
% (shown here in the data_unc_p1) variable.
% Several nuclei are merged together and a large portion of the
% bottom-left is surrounded by uncertainty where Ilastik could not
% determine exactly what was background and what was nuclei. 
% But the biggest issues are segmentation errors with the merged nuclei.
% evidenced by the image on variable t_p1 which also coincides with the
% region of biggest undertainty.

% improving mask from part1

data_seg_p1_mask = data_seg_p1 > 0.9;
imshow(data_seg_p1_mask);

% eroding
s = round(1.15*sqrt(mean(area))/pi);
nucmin = imerode(fusedMask, strel('disk',s));
imshow(nucmin, 'InitialMagnification', 'fit');

% get outside region
outside = ~imdilate(fusedMask, strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');

% basins for ws
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin);
shading flat;

L = watershed(basin);
rgb = label2rgb(L,'jet',[.5 .5 .5]);
imshow(rgb);


% combining
newNuclearMask = L > 11 | (data_seg_p1_mask - fusedMask);
imshow(newNuclearMask, 'InitialMagnification', 'fit');

% comparing Ilastik segmentation from part 1 with improved mask on matlab
imshowpair(newNuclearMask, data_seg_p1);

% Even though using a watershed function in matlab was able to split some
% of the large merged areas in the figure we can still visualize several
% problems. The area on the bottom left corner is still not well divided
% and several nucleis are still merged together. A small area with smaller
% nuclei had all nuclei dissapear due to the effect of erosion. Since the
% uncertainty in that area was too big, increasing the erosion would result
% in an increase in the number of nuclei that would dissapear. 
% Matlab improvement did a better job in central areas of the figure were
% uncertainty from Ilastik was smaller and the segmentation from Ilastik
% was already better


% Part 3. Redo part 1 but now be more aggresive in defining the background.
% Try your best to use ilastik to separate cells that are touching. Output
% the resulting mask into the repository. What is the problem now?

% part 3 and part 4 will be answered together below

% Part 4. Read your mask from Part 3 into MATLAB and try to improve
% it as best you can.

data_seg_p2 = h5read('Segmentation_Nuclei_part2.h5', '/exported_data');
data_seg_p2 = squeeze(data_seg_p2);
imshow(data_seg_p2, []);

% uncertainty data from the Ilastik analysis from part 1
data_unc_p2 = h5read('Uncertainty_part2.h5', '/exported_data');
data_unc_p2 = squeeze(data_unc_p2);
imshow(data_unc_p2, []);

% As we can see from the image above. If we are more aggressive to use
% ilastik to separate cells that are touching, the software does a much
% better job at separating the nuclei as we can clearly see in data_seg_p2, specially comparing the bottom left corner of part1 and part2.
% However, some corners of the nuclei are missing because now all the area
% of uncertainty is on the edges of the nuclei. Also, some nuclei present
% holes in their composition.

% improving mask from part2

data_seg_p2_mask = data_seg_p2 > 0.9;
imshow(data_seg_p2_mask);

data_seg_p2_dil = imopen(data_seg_p2_mask, strel('disk',3));
imshow(data_seg_p2_dil);

% a simple erosion followed by a dilation is enough to improve on the mask
% and fill the holes that exist in the data.

imshowpair(data_seg_p2_mask, data_seg_p2_dil);

%% Problem 2. Segmentation problems.

% The folder segmentationData has 4 very different images. Use
% whatever tools you like to try to segement the objects the best you can. Put your code and
% output masks in the repository. If you use Ilastik as an intermediate
% step put the output from ilastik in your repository as well as an .h5
% file. Put code here that will allow for viewing of each image together
% with your final segmentation. 

% Bacteria.tif

addpath('segmentationData/');

bacteria_orig = imread('bacteria.tif');
imshow(bacteria_orig);

bac = h5read('Segmentation_bacteria.h5', '/exported_data');
bac = squeeze(bac);
imshow(bac, []);

bac_d = ~imdilate(bac, strel('disk',3));
imshow(bac_d, []);

% segmentation of bacteria.tif was performed in Ilastik. Ilastik did a
% great job in segmenting the file. imdilate was used to fill tiny holes in
% the resulting segmented image.

%cellPhaseContrast.tif

cell_orig = imread('cellPhaseContrast.png');
imshow(cell_orig);

cell_seg = h5read('Segmentation_cell.h5', '/exported_data');
cell_seg = squeeze(cell_seg);
cell_seg = imrotate(cell_seg,270);
imshow(cell_seg, []);

%watershed

cell_seg_mask = cell_seg > 0.2;
imshow(cell_seg_mask);

CC = bwconncomp(cell_seg_mask);
stats = regionprops(CC, 'Area');
area = [stats.Area];
fusedCandidates = area > mean(area) + std(area);
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1, sublist{:});
fusedMask = false(size(cell_seg_mask));
fusedMask(sublist) = 1;
imshow(fusedMask, 'InitialMagnification', 'fit');

% eroding
s = round(1.2*sqrt(mean(area))/pi);
nucmin = imerode(fusedMask, strel('disk',s));
imshow(nucmin, 'InitialMagnification', 'fit');

% get outside region
outside = ~imdilate(fusedMask, strel('disk',1));
imshow(outside, 'InitialMagnification', 'fit');

% basins for ws
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);
pcolor(basin);
shading flat;

L = watershed(basin);
rgb = label2rgb(L,'jet',[.5 .5 .5]);
imshow(rgb);

% combining
newNuclearMask = L > 1 | (cell_seg_mask - fusedMask);
imshow(newNuclearMask, 'InitialMagnification', 'fit');

% ilastik was used to perform the initial segmentation. Ilastik was good
% for detecting the initial nuclei, however it was not good at separating
% the fused nuclei due to their similarity color with the background in
% phase contrast therefore watersheding was used to improve the separation
% of the fused nuclei.

% worms.tif

worms_orig = imread('worms.tif');
imshow(worms_orig);

worms_seg = h5read('Segmentation_worms.h5', '/exported_data');
worms_seg = squeeze(worms_seg);
worms_seg = imrotate(worms_seg,270);
worms_seg = flipdim(worms_seg ,2);
imshow(worms_seg, []);

worms_open = ~imopen(worms_seg, strel('disk',4));
imshow(worms_open);

% ilastik was used to perform the initial segmentation and background removal. Ilastik
% segmentation was already good, just a simple erosion followed by dilation
% was enough to improve the segmentation

%yeast.tif

yeast_orig = imread('yeast.tif');
imshow(yeast_orig);

BW = im2bw(yeast_orig, graythresh(yeast_orig));
imshow(BW);

% a simple conversion to black and white creates a good mask of the image
% of the yeasts

