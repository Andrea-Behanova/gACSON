function [J,inxCollection,flag] = util_rgn_superVox(I,ini_inx,t_dist,area,fg)

[r,c,h] = size(I);
J = false(r,c,h); J(ini_inx) = 1;
flag = 0;

if fg == 1
    ngb = [r;-r;1;-1;(r*c);-(r*c)];
else
    ngb1 = [r-1; r; r+1; -r+1; -r; -r-1; 1; -1;];
    ngb2 = ngb1 - (r*c);
    ngb3 = ngb1 + (r*c);
    ngb = [ngb2;ngb1;ngb3];
end
nngb = length(ngb);
L_ini_inx = length(ini_inx);
rgn_mean = sum(I(ini_inx))/L_ini_inx;
inx = zeros(L_ini_inx,nngb);
inxCollection = zeros(area(2),1);
L_rgn = 1;
inxCollection(1:L_ini_inx) = ini_inx;

while L_ini_inx ~= 0
    rgn_mean_old = rgn_mean;
    for j = 1:nngb
        inx(:,j) = ini_inx + ngb(j);
    end
    inx = unique(inx);
    dist = abs(I(inx) - rgn_mean);
    ini_inx = inx((dist<t_dist) & J(inx)==0);
    J(ini_inx) = 1;
    L_ini_inx = length(ini_inx);
    old_intn = L_rgn*rgn_mean_old;
    inxCollection(L_rgn+1:L_rgn+L_ini_inx) =  ini_inx;
    L_rgn = L_rgn + L_ini_inx;
    if L_rgn>area(2); flag = 2; break; end;
    rgn_mean = (old_intn+sum(I(ini_inx))) / L_rgn;
    inx = zeros(L_ini_inx,nngb);
end
if L_rgn<area(1); flag = 1; end;
inxCollection = inxCollection(inxCollection~=0);