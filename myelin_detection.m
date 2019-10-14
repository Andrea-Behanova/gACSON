function  myelin_detection(pic,handles)

slmin = 1;
slmax = 256;
thr = 70;

f = pic;
f(f<=thr)=true;
f(f>thr)=false;
f = logical(f);
f = bwareaopen(f, 500);
fig = figure;
imshow(f)
%set(fig, 'Position', [1000 600 700 700]);
S = [];
S.thr = thr;


h = uicontrol('Style', 'slider', 'Callback', {@sliderCallback,handles}, 'Min',slmin,'Max',slmax,...
                'SliderStep',[1 1]./(slmax-slmin),'Value',70,...
                'Position',[20 20 200 20]);
% h.image = pic;
%     f = pic;
%     f(f<thr)=0;
%     f(f>thr)=1;
%     imshow(f)
function thr = sliderCallback(hObject, evt, handles)
    thr = get(hObject,'Value');
%     a = get(gca,'children');
%     pic = a.CData;
    f = pic;
    f(f<=thr)=true;
    f(f>thr)=false;
    f = logical(f);
    f = bwareaopen(f, 200);
    set(get(gca,'children'),'CData',f)
    assignin('base','threshold',thr)
end

S.pb = uicontrol('style','push',...
                 'units','pix',...
                 'position',[400 20 200 40],...
                 'fontsize',14,...
                 'string','Selection finished',...
                 'callback',{@pb_call,S,handles});
             
 function pb_call(varargin)
% S = varargin{3};  % Get the structure.

close(fig)
 end


end