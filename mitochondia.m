clear all
%% detection of mitochondrias
% filename = 'lbl_axon_160520_Contra_1.mat';
% load(filename)
% label = lbl;
% refined_seg = label; 
% 
% mAxon_lbl = unique(refined_seg);
% mAxon_lbl = nonzeros(mAxon_lbl);
% 
% load('fil_im_160520_Contra_1.mat')
% I = filt_im;
% sv = 11;
% opts.supervoxelsize= [sv sv sv];
% opts.spacing = [1 1 3.6];
% opts.compactness = 23;
% opts.numIter = 5;
% I = (round(mat2gray(I)*255));
% S = SLICSupervoxelsMex(uint32(I),opts);
% s = size(I);
% S = reshape(S,s);
% SV = S+1;
% 
% 
% 
% fibers = util_mitochondria_supervoxels(refined_seg,mAxon_lbl,SV);
% save(['lbl_mito_' filename(10:end)],'fibers','-v7.3')
% 
% %label = imadd(refined_seg,fibers);
% %save(filename,'label','-v7.3')

%% connect axon with mitochondrias
filename = 'lbl_axon_160506_Ipsi_1.mat';
load(filename)
axon = label;

load('lbl_mito_160506_Ipsi_1.mat');
mito = lbl;

axon(logical(mito)) = mito(logical(mito));

save(filename,'axon','-v7.3')