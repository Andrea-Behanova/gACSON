function f = features(hObject, eventdata,handles,im,pos)
%gaussian smoothing
sigma_gs = 2;
GS = imgaussfilt(im,sigma_gs);

%laplacian of gaussian
sigma_log = 0.5;
w=fspecial('log',[3 3],sigma_log); 
LoG=imfilter(im,w,'replicate'); 
% imshow(LoG);

%gaussian gradient magnitude
sigma_ggm = 2;
GS2 = imgaussfilt(im,sigma_ggm);
[Gx,Gy] = imgradientxy(GS2);
[GGm,Gdir] = imgradient(Gx,Gy);
% figure
% imshow(GGm);

%difference of gaussian
sigma_dog = 5;
gauss = fspecial('gaussian', [10 10], sigma_dog);
blur1 = imfilter(im, gauss, 'replicate');
DoG = im - blur1;
% figure


linPos1 = sub2ind(size(im), pos(:,2), pos(:,1));
f1 = GS(linPos1);
f2 = LoG(linPos1);
f3 = GGm(linPos1);
f4 = DoG(linPos1);

f = [f1 f2 f3 f4];