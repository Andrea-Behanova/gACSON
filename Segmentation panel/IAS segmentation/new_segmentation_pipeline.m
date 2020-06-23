function new_segmentation_pipeline(opt,hObject,handles)
tic
g = waitbar(0,'Watershed centroids');


r_im = opt.read_im;
% s_address = opt.save_address;
r_mask = opt.read_mask;
myelin_inx = opt.myelin_inx;
similarity_th = opt.similarity_th;
min_max_myelin_volume = opt.min_max_myelin_volume;
min_max_lbl_volume = opt.min_max_lbl_volume;
small_volume_th = opt.closing_small_volume_threshold;
bg = opt.background;
fg = opt.foreground;
%% Loads

% Load stack of raw images after alignmet.
% t = load(r_im);

indx = strcmp(handles.listFiles.String, r_im);
indx = indx(2:end);
handles.indx = indx;
image = handles.image(indx);
raw_im = image{1,1};
% raw_im = raw_im(:,:,1:300);

% fields = fieldnames(t); raw_im = t.(fields{1});
% raw_im = squeeze(raw_im);

% Load mask
% AFTER alignment, non-tissue voxels are filled with intensity=255.
% Mask is equal to one for the non-tissue voxels.
if isempty(r_mask)
    mask = [];
else
    t = load(r_mask);
    fields = fieldnames(t); mask = t.(fields{1});
    mask = mask(:,:,1:300);
end
clear t fields

waitbar(0.1,g);
%% BM4D filtering

if opt.filt_im == 0
    [filt_im,~] = bm4d(raw_im,'Gauss',0, 'lc');
%     save(strcat('C:\Users\andrb\Desktop\imshow3D\GUIsegmentation\','filt_im_160323'),'filt_im','-v7.3')
    clear raw_im
else
    filt_im = raw_im;
    clear raw_im  
end

%save('fil_im_151113_ipsi_1','filt_im','-v7.3')

filtered_image = filt_im;
filtered_image2 = filt_im;
waitbar(0.2,g);
%% Canny edge detection
% STD Gaussian filter: sqrt(2).
% Weak and strong edges:0.25 and 0.6 times the maximum gradient magnitude.

filt_im = double(filt_im);
filt_im(mask) = 255;
bound_im = false(size(filt_im));
for i = 1:size(filt_im,3)
    bound_im(:,:,i) = edge(filt_im(:,:,i),'canny');
end
bound_im(mask) = true;
%save(strcat(s_address,'mat_edge'),'bound_im','-v7.3')

bound_im = imdilate(bound_im,true(3));

waitbar(0.3,g);
%% Myelin detection

if myelin_inx == 0
    thr = evalin('base', 'threshold');
    f = filt_im;
    f(f<=thr)=true;
    f(f>thr)=false;
    myelin = logical(f);
    myelin = bwareaopen(myelin, 2000);
    
elseif myelin_inx == 1
    myelin = handles.ML.label;
    
else
    filt_im = padarray(filt_im/255, [1 1 1], 1);
    bound_im = padarray(bound_im, [1 1 1], 1);
    [r1,c1,h1] = size(filt_im);

    myelin_inx = sub2ind([r1,c1,h1],myelin_inx(:,1),myelin_inx(:,2),myelin_inx(:,3));
    myelin_rgn = false(r1,c1,h1);

    if isinf(min_max_myelin_volume(2))
        min_max_myelin_volume(2) = numel(filt_im);
    end

    [myelin_rgn_inx,~] = util_regionGrowing(filt_im,myelin_inx,similarity_th,min_max_myelin_volume,bound_im,bg);
    myelin_rgn(myelin_rgn_inx) = true;

    myelin_rgn = imclose(myelin_rgn,true(3));

    myelin = myelin_rgn;
    myelin = myelin(2:end-1,2:end-1,2:end-1);
end

disp('Myelin detection finished')

gg = figure;
imshow(myelin(:,:,1))
answer = questdlg('Are you satisfyied with the result? Would you like to continue?', ...
	'Myelin detection', ...
	'Yes, continue','No, return to the programm','Yes, continue');
% Handle response
switch answer
    case 'Yes, continue'
        close(gg)
    case 'No, return to the programm'
        close(gg)
        return
        close(g)
end

waitbar(0.4,g);

%% Myelin supervoxels
I = filtered_image;

sv = 11;
opts.supervoxelsize= [sv sv sv];
opts.spacing = [1 1 3.6];
opts.compactness = 23;
opts.numIter = 5;
I = (round(mat2gray(I)*255));

S = SLICSupervoxelsMex(uint32(I),opts);



s = size(I);
S = reshape(S,s);
S = S+1;

refinedMyelin = util_input_refine_axon(filtered_image, uint8(myelin),  S, [1 1 50/15], 'm');

myelin_rgn = padarray(refinedMyelin, [1 1 1], 1);
disp('Supervoxels myelin finished')

waitbar(0.5,g);
%% watershed in 2D
seed = [];

f = waitbar(0,'Seeds reduction');
bound_im = myelin_rgn;
for i = 1:size(myelin_rgn,3)
    pic = myelin_rgn(:,:,i);

    C = pic;
    C = ~bwareaopen(~C, 100);
    D = -bwdist(C);

    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    s = regionprops(Ld2,'Centroid');

    d1 = imdilate(Ld2,true(17))-Ld2;
    pic(logical(d1))=1;
    bound_im(:,:,i) = pic;

    centroids = floor(cat(1, s.Centroid));

    lin_inx = sub2ind(size(myelin_rgn),centroids(:,2),centroids(:,1), repmat(i,size(centroids,1),1));
    % h(h(:,2)>size(Ld2,2))=size(Ld2,2);
    % h(h(:,1)>size(Ld2,1))=size(Ld2,1);

    seed = [seed;lin_inx];
    frac = i/(size(myelin_rgn,3)/100)/100;
    waitbar(frac,f)
end

waitbar(0.6,g);

disp('Seeds reduction finished')
bound_im = bound_im(2:end-1,2:end-1,2:end-1);

min_max_lbl_volume = [100 2000000];

% bound_im = padarray(bound_im, [1 1 1], 1);
bw_close = bwareaopen(~bound_im,small_volume_th);
small_volumes = ~(bw_close|bound_im);
bound_im(small_volumes) = true;

filtered_image = padarray(filtered_image/255, [1 1 1], 1);
bound_im = padarray(bound_im, [1 1 1], 1);
[r1,c1,h1] = size(filtered_image);

rps = seed(2:end)';

label = zeros(r1,c1,h1);
lbl = 4;
h = 0;

for seed1 = rps
    if bound_im(seed1)==0
%           seed1 = sub2ind(size(myelin_rgn),50,411, 185);
        [rgn_inx,flag] = util_regionGrowing(filtered_image,seed1,similarity_th,min_max_lbl_volume,bound_im,fg);
        if  flag==0            
            label(rgn_inx) = lbl;
            bound_im(rgn_inx) = true;
            lbl = lbl+1;
        end
    end
    h = h + 1/length(rps);
    waitbar(h,f,'Seeds growing')
end

label = label(2:end-1,2:end-1,2:end-1);
disp('Seeds growing finished')
waitbar(0.7,g);
%% SUPERVOXELS + removing oversegmented intraaxons
I = filtered_image2;
% [r, c, sz] = size(I);

sv = 11;
s = 1;
opts.supervoxelsize= [sv sv sv];
opts.spacing = [1 1 3.6];
opts.compactness = 23;
opts.numIter = 5;
I = (round(mat2gray(I)*255));

S = SLICSupervoxelsMex(uint32(I),opts);
disp('Supervoxels finished')


s = size(I);
S = reshape(S,s);
S = S+1;

refinedIntraAxon = util_input_refine_axon(filtered_image2, label,  S, [1 1 50/15], 'ias');

disp('Removing oversegmented axons finished')
waitbar(0.8,g);
%% removing unmyelinated axons
myelin_rgn = imdilate(refinedMyelin,true(8));
label = refinedIntraAxon;
for i = 1:max(label(:))
    %1:max(label(:))
    a = false(size(label,1),size(label,2),size(label,3));
    a(label==i) = 1;
    if sum(a(:))>500
        BW = imdilate(a,true(9));
        BW(label>0)=false;
%         BW = edge3(a,'approxcanny',0.5); 
        if sum(myelin_rgn(BW))<0.4*sum(BW(:))
            label(label==i) = 0;
        end
    end
    frac = i/(max(label(:))/100)/100;
    waitbar(frac,f,'Removing unmyelinated axons');
end
disp('Removing unmyelinated axons finished')
save(['lbl_axon' r_im(7:end)],'label','-v7.3')
waitbar(0.9,g);
%% watershed myelin
lbl = false(size(label,1),size(label,2),size(label,3));
lbl(label>0)=true;
lbl = imfill(lbl,'holes'); 
lbl = ~bwareaopen(~lbl, 500);

D = bwdist(lbl);

Ld2 = watershed(D);

refinedMyelin(Ld2==0)=0;

myelin_lbl = immultiply(logical(refinedMyelin),Ld2);
disp('Watershed myelin finished')
waitbar(0.9,g);
save(['lbl_myelin' r_im(7:end)],'myelin_lbl','-v7.3')
toc
waitbar(1,g);
close(f)
close(g)
