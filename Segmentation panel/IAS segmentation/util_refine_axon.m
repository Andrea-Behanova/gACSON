function [refinedAxon,eqSV] = util_refine_axon(cropIm,cropAxon,cropSV)

cropSVIndx = label2idx(cropSV);
coveredSuperVoxel = double(cropSV).*cropAxon;
coveredSuperVoxelIndx = label2idx(coveredSuperVoxel);
eqSV = false(1,length(cropSVIndx));
flag = 0;
newAxon = zeros(size(cropIm));
for i = 1:length(coveredSuperVoxelIndx)
    if ~isempty(coveredSuperVoxelIndx{i})
        a = length(coveredSuperVoxelIndx{i});
        b = length(cropSVIndx{i});
        r = a/b;
        if r>0.3
            eqSV(i) = true;
            newAxon(cropSVIndx{i}) = mean(cropIm(cropSVIndx{i}));
            flag = 1;
        end
    end
end

if flag==1
    exe_im1 = padarray(newAxon/255, [1 1 1], inf);
    exe_im2 = padarray(cropAxon, [1 1 1]);
    axon_inx = find(exe_im2&exe_im1);
    [refinedAxon,~,~] = util_rgn_superVox(exe_im1,axon_inx,0.1,[0,1000000],1);
    refinedAxon = refinedAxon(2:end-1,2:end-1,2:end-1);
else
    refinedAxon = cropAxon;
end
refinedAxon = refinedAxon | cropAxon;


