function util_HM_SBEM_segmentation_pipeline(opt,hObject,handles)
g = waitbar(0,'Regional maxima');

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
    [filt_im,~] = bm4d(raw_im,'Gauss',0, handles.opt.filt_prof);
%     save(strcat('C:\Users\andrb\Desktop\imshow3D\GUIsegmentation\','filt_im_160323'),'filt_im','-v7.3')
    clear raw_im
else
    filt_im = raw_im;
    clear raw_im  
end

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
    
    bound_im(myelin) = true;
    
    filt_im = padarray(filt_im/255, [1 1 1], 1);
    bound_im = padarray(bound_im, [1 1 1], 1);
    myelin_rgn = padarray(myelin, [1 1 1], 1);
    [r1,c1,h1] = size(filt_im);
    label = double(myelin_rgn);
    
elseif myelin_inx == 1
    myelin = handles.ML.label;
    
    bound_im(myelin) = true;
    
    filt_im = padarray(filt_im/255, [1 1 1], 1);
    bound_im = padarray(bound_im, [1 1 1], 1);
    myelin_rgn = padarray(myelin, [1 1 1], 1);
    [r1,c1,h1] = size(filt_im);
    label = double(myelin_rgn);
    
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
    
    bound_im(myelin_rgn) = true;
    label = double(myelin_rgn);

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
end

waitbar(0.4,g);
%% Closing small volumes

bw_close = bwareaopen(~bound_im,small_volume_th);
small_volumes = ~(bw_close|bound_im);
bound_im(small_volumes) = true;

small_volumes = small_volumes(2:end-1,2:end-1,2:end-1);
%save(strcat(s_address,'mat_small_volumes'),'small_volumes','-v7.3')

clear small_volumes bw_close

waitbar(0.6,g);
%% Finding the location of seeds

rgnmaxima = zeros(r1,c1,h1);
for i = 2:h1-1
    D = bwdist(bound_im(:,:,i));
    rgnmaxima(:,:,i) = imregionalmax(D,4);
end
rps = find(rgnmaxima);
rps = rps'; 
clear rgnmaxima

waitbar(0.8,g);
%% Growing a volume for each seed
% lbl starts from an arbitrary value, here 4; in case we wanted to hold
% all small volumes, edges ... in the "label" matrix.     


lbl = 4;
tic
f = waitbar(0,'Seeds growing');

for seed = rps
    h = seed/rps(end);
    
    guidata(hObject, handles);

    if bound_im(seed)==0
        [rgn_inx,flag] = util_regionGrowing(filt_im,seed,similarity_th,min_max_lbl_volume,bound_im,fg);
        if  flag==0            
            label(rgn_inx) = lbl;
            bound_im(rgn_inx) = true;
            lbl = lbl+1;
        end
    end
    drawnow limitrate
    
    if seed == rps(1000)
        timer = toc;
        timer = timer/1000;
        app_time = round(timer*(length(rps))/3600,2);
        if app_time < 1
            app_time = [num2str(round(app_time*60)) ' minutes'];
        else
            app_time = [num2str(app_time) ' hours'];
        end
        answer = questdlg(['Approximation time is: ' app_time ' Would you like to continue?'], ...
        'Time approximation', ...
        'Yes','No','Yes');
        % Handle response
        switch answer
            case 'Yes'
                continue
            case 'No'
%                 set(handles.perc_string, 'String', [0 ' %'])
%                 set(handles.percentage,'Position', [0.19 0.01 0 0.03]);
                guidata(hObject, handles);
                return
        end
    end
    waitbar(h,f);
end
    drawnow
waitbar(1,g);
label = label(2:end-1,2:end-1,2:end-1);
[filename, filepath] = uiputfile('*.mat', 'Save the project file:');
FileName = fullfile(filepath, filename);
save(FileName, 'label', '-v7.3');

close(f)
close(g)
