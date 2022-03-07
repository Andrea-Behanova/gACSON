function [label] = Choosed_seed_segmenttaion(opt,hObject,handles)

[y,x] = getpts; y = round(y(:,1)); x = round(x(:,1));
z = (handles.S+1)*ones(1,length(x))';
choice = [x+1, y+1, z];
g = waitbar(0,'Manual selection');

r_im = opt.read_im;
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
end
clear t fields

waitbar(0.1,g);
%% BM4D filtering

if opt.filt_im == 0
    [filt_im,~] = bm4d(raw_im);
    save(strcat(s_address,'filt_im'),'filt_im','-v7.3')
    clear raw_im
else
    filt_im = raw_im;
    clear raw_im  
end


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
    
    bound_im(myelin) = true;
    
    filt_im = padarray(filt_im/255, [1 1 1], 1);
    bound_im = padarray(bound_im, [1 1 1], 1);
    [r1,c1,h1] = size(filt_im);
    
elseif myelin_inx == 1
    myelin = handles.ML.label;
    
    bound_im(myelin) = true;
    
    filt_im = padarray(filt_im/255, [1 1 1], 1);
    bound_im = padarray(bound_im, [1 1 1], 1);
    [r1,c1,h1] = size(filt_im);
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
    % label = double(myelin_rgn);

    bound_im(myelin_rgn) = true;

    myelin_rgn = myelin_rgn(2:end-1,2:end-1,2:end-1);
    %save(strcat(s_address,'mat_myelin_rgn'),'myelin_rgn','-v7.3')

    clear myelin_rgn
end

waitbar(0.4,g);
%% Closing small volumes

bw_close = bwareaopen(~bound_im,small_volume_th);
small_volumes = ~(bw_close|bound_im);
bound_im(small_volumes) = true;

small_volumes = small_volumes(2:end-1,2:end-1,2:end-1);
%save(strcat(s_address,'mat_small_volumes'),'small_volumes','-v7.3')

clear small_volumes bw_close


%% Finding the location of seeds

rgnmaxima = zeros(r1,c1,h1);
for i = 1:size(choice,1)
    rgnmaxima(choice(i,1),choice(i,2),choice(i,3)) = 1;
end
rps = find(rgnmaxima);
rps = rps'; 
clear rgnmaxima

waitbar(0.6,g);
%% Growing a volume for each seed
% lbl starts from an arbitrary value, here 4; in case we wanted to hold
% all small volumes, edges ... in the "label" matrix.     

label = zeros(r1,c1,h1);
lbl = 1;
for seed = rps
    if bound_im(seed)==0
        [rgn_inx,flag] = util_regionGrowing(filt_im,seed,similarity_th,min_max_lbl_volume,bound_im,fg);
        if  flag==0            
            label(rgn_inx) = lbl;
            bound_im(rgn_inx) = true;
            lbl = lbl+1;
        end
    end
    
end

label = label(2:end-1,2:end-1,2:end-1);
waitbar(1,g);
close(g)