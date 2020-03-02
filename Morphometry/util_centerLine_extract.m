function ordered_CentLine = util_centerLine_extract(skel,sz)

bwSkel = false(sz);
for i = 1:length(skel)
    L = ceil(skel{i});
    for j = 1:size(L,1)
        linearInd = sub2ind(sz, L(j,1), L(j,2), L(j,3));
        bwSkel(linearInd) = 1;
    end
end

[lbl,n] = bwlabeln(bwSkel,26); 
if n > 1
    stats = regionprops(lbl,'Area');
    [~,inxMain] = max([stats.Area]);
    bwSkel = lbl == inxMain;
end


[bwCentLine,sp] = util_skeleton_shortest_path(bwSkel);
lsp = sub2ind(sz, sp(1), sp(2), sp(3));
g = bwdistgeodesic(bwSkel,lsp);

branch = util_bw2skel_equv(skel,bwCentLine);

accurateCentLine = cell(1);
j = 1;
for i = 1:length(branch)
    if branch(i)
        accurateCentLine(j) = skel(i);       
        L = ceil(skel{i});
        ls = sub2ind(sz, L(1,1), L(1,2), L(1,3)); le = sub2ind(sz, L(end,1), L(end,2), L(end,3));        
        se_dist(j,:) = [g(ls),g(le)];
        j = j+1;
    end
end

j = 1;
ordered_CentLine = cell(1);
b_inx = find(branch);
c_se_dist = se_dist;
d = unique(se_dist);
for i = 1:length(d)
    [r,c] = find(c_se_dist == d(i));
    if length(r)>1
        dis = zeros(length(r),1);
        for jj = 1:length(r)
            dis(jj) = abs(se_dist(r(jj),1)-se_dist(r(jj),2));
        end
        [~,rinx] = max(dis);
        tt = r;
        tt(rinx) = [];
        c_se_dist(tt,:) = nan;
        r = r(rinx);
    end
    c_se_dist(r,:) = nan;
    if c==1
        ordered_CentLine(j) = skel(b_inx(r));
    elseif c==2
        ordered_CentLine(j) = {flipud(skel{b_inx(r)})};
    end
    j = j+1; 
end


