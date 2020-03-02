function fibers = util_mitochondria_supervoxels(refined_seg,mAxon_lbl,SV,filename)

sz = size(refined_seg);

s3D = regionprops(refined_seg, 'BoundingBox', 'PixelIdxList');
pixelIndexList = {s3D.PixelIdxList};

fibers = zeros(sz);
for N = mAxon_lbl'
    lbl = false(sz);
    lbl(pixelIndexList{N}) = true;
    bb = s3D(N).BoundingBox;
    [cropAxon,imb] = util_extract_bounded_obj(lbl,bb,[10,10,10]);
    cropSV = util_extract_bounded_obj(SV,bb,[10 10 10]);
    
    p = 45; padCropAxon = padarray(cropAxon,[p,p,p]);
    closed_cropAxon = false(size(cropAxon));
    searchRGN = bwdist(padCropAxon); searchRGN = searchRGN<40;
    dSearchRGN = bwdist(~searchRGN); dSearchRGN = dSearchRGN>40;
    dSearchRGN = dSearchRGN(p+1:end-p,p+1:end-p,p+1:end-p);     
    closed_cropAxon(dSearchRGN) = 1;
    mito = closed_cropAxon & ~cropAxon;
    refined_mito_inx = util_sv_equal(mito,cropSV);
    
    [r,c,h] = ind2sub(size(mito),refined_mito_inx);
    nr = r+imb(1)-1; nc = c+imb(3)-1; nh = h+imb(5)-1;
    inx=sub2ind(sz, nr, nc, nh);
    fibers(inx) = N;
end
%save(strcat(filename,'fibers'),'fibers','-v7.3')


function refined_mito_inx = util_sv_equal(mito,cropSV)

svlbl = nonzeros(unique(cropSV(mito)));
pixelIndexList_SV = label2idx(cropSV);
sv_volume = cellfun(@(x) length(x), pixelIndexList_SV(svlbl));
mito_volume_of_sv = cellfun(@(x) sum(mito(x)), pixelIndexList_SV(svlbl));
mito_sv = svlbl((mito_volume_of_sv./sv_volume) > 0.6);
refined_mito_inx = cell2mat({pixelIndexList_SV{mito_sv}}');



