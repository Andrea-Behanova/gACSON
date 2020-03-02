function [centLine,s] = util_skeleton_shortest_path(bwSkel)

bwSkel = padarray(bwSkel, [1 1 1]);
D = bwdist(~bwSkel);

[~,indmaxD] = max(D(:));

% bwSkel(indmaxD);

geoD = bwdistgeodesic(bwSkel,indmaxD);
[~,sp] = max(geoD(:));

[s(1),s(2),s(3)] = ind2sub(size(bwSkel),sp); s = s-1;

bwdg = bwdistgeodesic(bwSkel,sp);

startPoint = find(bwdg==max(bwdg(:)));
startPoint = startPoint(1);
bwdg(startPoint) = NaN;

centLine = false(size(bwSkel));
centLine(startPoint) = true;

[r,c,~] = size(bwSkel);
ngb1 = [r-1; r; r+1; -r+1; -r; -r-1; 1; -1];
ngb2 = ngb1 - (r*c);
ngb3 = ngb1 + (r*c);
ngb4 = [-r*c, r*c];
neighborOffsets = [ngb2;ngb1;ngb3];
neighborOffsets = [neighborOffsets(:);ngb4(:)];

neighbors = bsxfun(@plus,startPoint,neighborOffsets);
[sortedNeighbors,inds] = sort(bwdg(neighbors));
activePixel = neighbors(inds(1));

% p = patch(isosurface(bwSkel,0.5)); p.FaceColor = 'red'; p.EdgeColor = 'none'; p.FaceAlpha = 0.3;

minRep(1) = max(bwdg(:));
i = 1;
while minRep(i)~=0
    
    centLine(activePixel)= true;
    bwdg(activePixel) = NaN;
    
%     [pt(1),pt(2),pt(3)] = ind2sub(size(bwSkel),activePixel);   
%     hold on
%     plot3(pt(2),pt(1),pt(3),'b*');
%     axis equal

    neighbors = bsxfun(@plus,activePixel,neighborOffsets);

%     max(bwdg(neighbors));
    activePixel = neighbors(bwdg(neighbors)==min(bwdg(neighbors)));
    i = i + 1;
    minRep(i) = min(bwdg(neighbors));
    if (~isempty(activePixel)) && (minRep(i) < minRep(i-1))
        activePixel = activePixel(1);
    elseif minRep(i-1) ~= 0
        activePixel = find(bwdg == minRep(i-1)-1);
        if ~isempty(activePixel)
            activePixel = activePixel(1);
            minRep(i) = minRep(i-1)-1;
        else
            break
        end
    end
end
centLine = centLine(2:end-1,2:end-1,2:end-1);




