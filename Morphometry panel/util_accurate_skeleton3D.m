function [skel,nBW] = util_accurate_skeleton3D(bw,res)

res = res./res(1);
[r,c,~] = size(bw);
r1 = ceil(1/res(3)*r); c1 = ceil(1/res(3)*c); 
rsBW = imresize(bw,[r1,c1],'box');

lbl = bwlabeln(rsBW,6);
stats = regionprops(lbl,'Area');
[~,inxMain] = max([stats.Area]);
if ~isempty(inxMain)
    nBW = lbl==inxMain;
    z = find(sum(sum(nBW,1),2));
    skel = zeros(numel(z),3);
    for jj = z'
        im_2d = nBW(:,:,jj);
        [I, J] = ind2sub([r1, c1],find(im_2d));
        skel(jj, :) = [mean(I), mean(J), jj];
    end
    skel = movmean(skel, 11);
%     skel = skeleton(nBW,false);
else
    skel = [];
    nBW = [];
end
