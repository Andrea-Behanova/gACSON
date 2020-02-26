function branch = util_bw2skel_equv(skel,bwCentLine)
branch = false(1,length(skel));
inxBWCentLine = find(bwCentLine);
for i = 1:length(skel)
    L = ceil(skel{i});
    linearL = sub2ind(size(bwCentLine),L(:,1),L(:,2),L(:,3));
    s = sum(ismember(inxBWCentLine,linearL))/length(linearL);
    if s>0.2
        branch(i) = true;
    end
end