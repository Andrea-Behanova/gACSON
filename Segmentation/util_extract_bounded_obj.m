function [bbw,bound] = util_extract_bounded_obj(bw,s3D,offset)

[r,c,h] = size(bw);
sR = ceil(s3D(2))-offset(1);
if sR<1
    sR = 1;
end
    
lR = sR+s3D(5)-1+2*offset(1);
if lR>r
    lR = r;
end

sC = ceil(s3D(1))-offset(2);
if sC<1
    sC = 1;
end


lC = sC+s3D(4)-1+2*offset(2);
if lC>c
    lC = c;
end

sH = ceil(s3D(3))-offset(3);
if sH<1
    sH = 1;
end

lH = sH+s3D(6)-1+2*offset(3);
if lH>h
    lH = h;
end

bbw = bw(sR:lR,sC:lC,sH:lH);
bound = [sR,lR,sC,lC,sH,lH];
    

