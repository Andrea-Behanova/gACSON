function refinedIntraAxon = util_input_refine_axon(im,aclbl,superVoxelSLIC,res,type)

pixelIndexList = label2idx(aclbl);


volume = zeros(1,length(pixelIndexList));
for i = 1:length(pixelIndexList)
    volume(i) = length(pixelIndexList{i})*prod(res);
end


ptn_intraAxon_inx = find(volume>1000);
intraAxon_inx = false(1,length(ptn_intraAxon_inx));
intraAxon = zeros(size(aclbl));
for i = 1:length(ptn_intraAxon_inx)
    lbl_int = im(pixelIndexList{ptn_intraAxon_inx(i)});
    mean_lbl = mean(lbl_int(:));
    if strcmp(type,'m')
        if mean_lbl<100 %myelin<100 intraaxonal>140
            intraAxon_inx(i) = 1;
            intraAxon(pixelIndexList{ptn_intraAxon_inx(i)}) = ptn_intraAxon_inx(i);
        end
    elseif strcmp(type,'ias')
        if mean_lbl>140 %myelin<100 intraaxonal>140
            intraAxon_inx(i) = 1;
            intraAxon(pixelIndexList{ptn_intraAxon_inx(i)}) = ptn_intraAxon_inx(i);
        end
    end
end


s3D = regionprops(intraAxon,'BoundingBox');

refinedIntraAxon = zeros(size(im));
for N = 1:length(s3D)
    bb = s3D(N).BoundingBox;
    if bb(4)>1 && bb(5)>1 && bb(6)>1
        obj = false(size(im));
        obj(pixelIndexList{N}) = true;
        r = length(nonzeros(refinedIntraAxon(pixelIndexList{N})))/length(pixelIndexList{N});
        [cropAxon,imb] = util_extract_bounded_obj(obj,bb,[10 10 3]);
        cropim = util_extract_bounded_obj(im,bb,[10,10,3]);
        cropSV = util_extract_bounded_obj(superVoxelSLIC,bb,[10,10,3]);
        refinedAxon = util_refine_axon(cropim,cropAxon,cropSV);
        obj(imb(1):imb(2),imb(3):imb(4),imb(5):imb(6)) = refinedAxon;
        if r<0.2  
            refinedIntraAxon(obj) = N;
        else
            mergingInx = median(nonzeros(refinedIntraAxon(pixelIndexList{N})));
            refinedIntraAxon(obj) = mergingInx;
        end
    end
end
% save (filename,'refinedIntraAxon','-v7.3')

