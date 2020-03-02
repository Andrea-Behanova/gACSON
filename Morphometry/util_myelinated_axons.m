function myelinated_axons = util_myelinated_axons(raw_im,refined_seg,SLIC_SV,volume_th,intensity_th,mRatio_th,save_add)

sz = size(refined_seg);

s3D = regionprops(refined_seg,'BoundingBox', 'PixelIdxList');
pixelIndexList = {s3D.PixelIdxList};

myelin = false(sz);
myelin(pixelIndexList{1}) = true;

massive_lbls = cellfun(@(x) length(x)>volume_th, pixelIndexList);
mean_lbls = cellfun(@(x) mean(raw_im(x))>intensity_th, pixelIndexList);

th = massive_lbls & mean_lbls;

myelinated_axons = zeros(sz);
l = length(massive_lbls);
for N = 1:l
    if th(N)
        lbl = false(sz);
        lbl(pixelIndexList{N}) = true;
        bb = s3D(N).BoundingBox;
        cropMyelin = util_extract_bounded_obj(myelin,bb,[10 10 10]);
        cropAxon = util_extract_bounded_obj(lbl,bb,[10,10,10]);
        cropSV = util_extract_bounded_obj(SLIC_SV,bb,[10,10,10]);
        r = util_myelin_extraction(cropMyelin,cropAxon,cropSV);
        if r > mRatio_th
            myelinated_axons(pixelIndexList{N}) = N;
        end
    end
end
save(strcat(save_add,'mat_myelinated_axons'),'myelinated_axons','-v7.3')


function r = util_myelin_extraction(cropMyelin,cropAxon,cropSV)

SV_inx = label2idx(cropSV);
dilate_5_axons = imdilate(cropAxon,true(5));
dilate_3_axons = imerode(dilate_5_axons,true(3));
boundary = dilate_5_axons & ~dilate_3_axons;
b_SV = nonzeros(unique(cropSV(boundary)));
boundary_SV = false(size(cropAxon));
for i = 1:length(b_SV)
    boundary_SV(SV_inx{b_SV(i)}) = true;
end
intsct = cropMyelin & boundary_SV;
r = sum(intsct(:))/sum(boundary_SV(:));






