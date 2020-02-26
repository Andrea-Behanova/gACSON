function AllS2D = util_planeEq_alpha05(bw,skel,nTangVec,IAS,resolut)
% IAS = 1 =intraaxonal space
% IAS = 2 =myelin
AllS2D = struct('Point',[],'Area',[],'MajorAxisLength',[],'MinorAxisLength',[],'Eccentricity',[],'EquivDiameter',[],'Thickness',[]);
r = 30; res = 0.25; k = 1;
for i = 1:length(skel)
    L = skel{i};
    tngV = nTangVec{i};
   
    [x,y] = ndgrid(-r:res:r, -r:res:r); z = zeros(size(x)); xyz = [x(:),y(:),z(:)];

    for j = 1:size(L,1)
        
        point = L(j,:);
        AllS2D(k).Point = point;
        
        normal = tngV(j,:);

        phi = pi - acos(normal*[0 0 1]'/(norm(normal)*norm([0 0 1])));
        phi_deg = phi * (180/pi);
        normal_proj = cross(normal,[0 0 1])/norm(cross(normal,[0 0 1]));

        [XYZnew, ~, ~] = AxelRot(xyz', phi_deg, normal_proj, [0 0 0]);

        xx = XYZnew(1,:); yy = XYZnew(2,:); zz = XYZnew(3,:);

        xx = reshape(xx,size(x))+ point(1); 
        yy = reshape(yy,size(x))+ point(2); 
        zz = reshape(zz,size(x))+ point(3);
        
        cross_section = interp3(double(bw),yy,xx,zz);
        
%         slice(double(bw),yy,xx,zz); hold on ; p = patch(isosurface(smooth3(bw,'gaussian',15),0.5)); 
%         p.FaceColor = 'red'; p.EdgeColor = 'none'; p.FaceAlpha = 0.3;
%         set(gca,'PlotBoxAspectRatio',[1 1 3.6])
        if k==40
            k
        end


% slice(double(bw),yy,xx,zz); hold on ; p = patch(isosurface(smooth3(bw,'gaussian',15),0.5));
% p.FaceColor = 'red'; p.EdgeColor = 'none';
% set(gca,'PlotBoxAspectRatio',[1 1 3.6])
% p = patch(isosurface(smooth3(nBW2,'gaussian',15),0.5)); %nBW2 is loaded from saved file
% p.FaceColor = 'black'; p.EdgeColor = 'none'; p.FaceAlpha = 0.3;
% grid off
% camlight; material shiny; lighting gouraud;

        bw_cross_section = cross_section >= 0.5;
        if IAS == 1
            bw_cross_section = uint8(bw_cross_section);
        elseif IAS == 2 %myelin
            BW2 = false(size(bw_cross_section));
            CC = bwconncomp(bw_cross_section);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            if ~isempty(idx)
                BW2(CC.PixelIdxList{idx}) = 1;
                bw_cross_section = BW2;

    %             J = imerode(BW2,true(3));
    %             BW2 = BW2 - J;
    %             BW2 = bwlabel(BW2);
    %             Distan = bwdist(BW2~=1 & BW2~=0);
    %             BW2 = BW2.*Distan;
    %             thickness = median(nonzeros(BW2(:)));

                D = bwdist(~BW2);
                BW = bwskel(BW2);

                BW2 = immultiply(D,BW)*2;
                thickness = nonzeros(BW2(:));
%                 thickness = median(nonzeros(thickness));
                AllS2D(k).Thickness = thickness/4*resolut(3)/10^3; 
            end
        end
        
        s2D = regionprops(bw_cross_section,'Area','MajorAxisLength','MinorAxisLength','Eccentricity','EquivDiameter');
%         s = util_draw_ellipse(bw_cross_section,bw_cross_section);
%         drawnow
        if length(s2D) > 1
            [~,inx] = max([s2D(:).Area]);
            
            %added
%             if negativ == 1
%                 s2D(inx) = [];
%                 [~,inx] = max([s2D(:).Area]);
%             end
            %
            
            s2D = s2D(inx);
        end
        if ~isempty(s2D)
            AllS2D(k).Area = s2D.Area/16*resolut(3)^2/10^6;
            AllS2D(k).MajorAxisLength = s2D.MajorAxisLength/4*resolut(3)/10^3;
            AllS2D(k).MinorAxisLength = s2D.MinorAxisLength/4*resolut(3)/10^3;
            AllS2D(k).Eccentricity = s2D.Eccentricity;
            AllS2D(k).EquivDiameter = s2D.EquivDiameter/4*resolut(3)/10^3;
        else
            s2D = struct('Point',point,'Area',0,'MajorAxisLength',0,'MinorAxisLength',0,'Eccentricity',0,'EquivDiameter',0,'Thickness',0);
            AllS2D(k) = s2D;
        end
        k = k+1;
    end
end


