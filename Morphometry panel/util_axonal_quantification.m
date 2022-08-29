function [allS2D, allS2D2] = util_axonal_quantification(final_axon,res,myelin)

sz = size(final_axon);
s3D = regionprops(final_axon,'BoundingBox','PixelIdxList');
numLabel = length(s3D);

allSkel = cell(numLabel,1);
allorderedCentLine = cell(numLabel,1);
allQuants = cell(numLabel,1);
allS2D = [];
allS2D2 = [];

%myelin
s3D2 = regionprops(myelin,'BoundingBox','PixelIdxList');
numLabel2 = length(s3D2);
for N = 1:numLabel
    pixIdxList = s3D(N).PixelIdxList;
    if length(pixIdxList)>4000
        bb = s3D(N).BoundingBox;
        axon = false(sz);
        axon(pixIdxList) = true;
        [cropAxon,~] = util_extract_bounded_obj(axon,bb,[9,9,3]);
        for jj = 1:size(cropAxon, 3)
            im_2d = cropAxon(:,:,jj);
            im_2d = imfill(im_2d, 'holes');
            cropAxon(:,:,jj) = im_2d;
        end
        cropAxon = imclose(cropAxon,true(5,5,5));

        [skel,nBW] = util_accurate_skeleton3D(cropAxon,res);            
        allSkel(N) = {skel};
        
        %myelin
        pixIdxList2 = s3D2(N).PixelIdxList;
        bb2 = s3D2(N).BoundingBox;
        myelin2 = false(sz);
        myelin2(pixIdxList2) = true;
        [cropMyelin,~] = util_extract_bounded_obj(myelin2,bb2,[9,9,3]);
        cropMyelin = imclose(cropMyelin,true(5,5,5));

        myelin2 = false(size(cropMyelin));
        CC2 = bwconncomp(cropMyelin);
        numPixels2 = cellfun(@numel,CC2.PixelIdxList);
        [~,idx2] = max(numPixels2);
        myelin2(CC2.PixelIdxList{idx2}) = 1;
        
        [skel2,nBW2] = util_accurate_skeleton3D(myelin2,res); 
        %
                        
        if ~isempty(skel) && ~isempty(skel2)
%             orderedCentLine = util_centerLine_extract(skel,size(nBW));
%             allorderedCentLine(N) = {orderedCentLine};
            orderedCentLine = {skel};
%             nTangVec = cell(length(orderedCentLine),1);
%             for i = 1:length(orderedCentLine)    
%                 L = orderedCentLine(i);
%                 if ~isempty(L)
            nTangVec = {util_normalTangentVector(skel,0)};
%                 end
%             end
%             BW = false(size(nBW));
%             CC = bwconncomp(nBW);
%             numPixels = cellfun(@numel,CC.PixelIdxList);
%             [biggest,idx] = max(numPixels);
%             BW(CC.PixelIdxList{idx}) = 1;
            
%             BW2 = false(size(nBW2));
%             CC2 = bwconncomp(nBW2);
%             numPixels2 = cellfun(@numel,CC2.PixelIdxList);
%             [biggest,idx2] = max(numPixels2);
%             BW2(CC2.PixelIdxList{idx2}) = 1;
            
            allS2D = util_planeEq_alpha05(nBW,orderedCentLine,nTangVec,1,res);
            allS2D2 = util_planeEq_alpha05(nBW2,orderedCentLine,nTangVec,2,res);
            allQuants(N) = {allS2D};
        end
    end
    
    if size(allS2D2,2)>41
        allS2D2 = allS2D2(1,21:end-20);
        allS2D = allS2D(1,21:end-20);
    end
        
    
%     MajorAxisIAS = cat(1, allS2D.MajorAxisLength);
%     ThicknessMyelin = cat(1, allS2D2.Thickness);
%     Thickness = cat(1, allS2D2(21:(length(MajorAxisMyelin)-20)).Thickness);
%     r = []; 
%     for i = 21:(length(MajorAxisIAS)-20)
%         r = [r allS2D(i).MajorAxisLength*ones(1,length(allS2D2(i).Thickness))]; 
%     end
%     gRatio = MajorAxisIAS./(MajorAxisIAS+2*ThicknessMyelin);
%     figure
%     subplot(131)
%     hist(gRatio(21:end-20),50)
%     title('G-ratio')
%     subplot(132)
%     hist((ThicknessMyelin(21:end-20)/4)*50, 50)
%     title('Thickness')
%     subplot(133)
%     hist((MajorAxisIAS(21:end-20)/4)*50, 50)
%     title('Inner diameter')
%     
%     figure
%     plot((ThicknessMyelin(21:end-20)/4)*50)
%     hold on
%     plot((MajorAxisIAS(21:end-20)/4)*50,'r')
%     legend('Thickness','Inner diameter')
%     hold off
end

% save(strcat(save_address,'_allSkel'),'allSkel');
% save(strcat(save_address,'_allorderedCentLine'),'allorderedCentLine');
% save(strcat(save_address,'_quants'),'allQuants');

