function [result] = Random_forest_seg(hObject, eventdata,handles,im,slice)

    
%train
train = handles.MLseg.train;

X = single(train(:,1:end-1));
Y = train(:,end);

TreeObject=TreeBagger(100,X,Y,'Method','classification','NVarToSample','all');

%test
pos1 = handles.MLseg.lbl0;
if ~isempty(pos1)
    idxPos = pos1(:,3)==slice;
    pos1 = pos1(idxPos,1:2);
    linPos1 = sub2ind(size(im), pos1(:,2), pos1(:,1));
end

pos2 = handles.MLseg.lbl1;
if ~isempty(pos2)
    idxPos = pos2(:,3)==slice;
    pos2 = pos2(idxPos,1:2);
    linPos2 = sub2ind(size(im), pos2(:,2), pos2(:,1));
end


pos3 = handles.MLseg.lbl2;
if ~isempty(pos3)
    idxPos = pos3(:,3)==slice;
    pos3 = pos3(idxPos,1:2);
    linPos3 = sub2ind(size(im), pos3(:,2), pos3(:,1));
end

pos4 = handles.MLseg.lbl3;
if ~isempty(pos4)
    idxPos = pos4(:,3)==slice;
    pos4 = pos4(idxPos,1:2);
    linPos4 = sub2ind(size(im), pos4(:,2), pos4(:,1));
end

pos5 = handles.MLseg.lbl4;
if ~isempty(pos5)
    idxPos = pos5(:,3)==slice;
    pos5 = pos5(idxPos,1:2);
    linPos5 = sub2ind(size(im), pos5(:,2), pos5(:,1));
end

mask = ones(size(im));
if ~isempty(pos1)
    mask(linPos1) = 0;
end

if ~isempty(pos2)
    mask(linPos2) = 0;
end

if ~isempty(pos3)
    mask(linPos3) = 0;
end

if ~isempty(pos4)
    mask(linPos4) = 0;
end

if ~isempty(pos5)
    mask(linPos5) = 0;
end

test_pos = find(mask);
[tPos(:,2),tPos(:,1)] = ind2sub(size(im),test_pos);

test = features(hObject, eventdata,handles,im,tPos);
test = single(test);

tic
L = size(test,1);
step = round(L/50);
seq = 1:step:L;
YFIT = cell(1,10);

parfor i = 1:length(seq)-1
    YFIT{i} = predict(TreeObject,test(seq(i):(seq(i+1)-1),:));
end

YFIT{length(seq)} = predict(TreeObject,test(seq(end):L,:));
toc

pred = [];
for i = 1:length(YFIT)
    pred = [pred;YFIT{1,i}];
end

result = zeros(size(im));
result(test_pos) = str2num(cell2mat(pred));

if ~isempty(pos1)
    result(linPos1) = 0;
end

if ~isempty(pos2)
    result(linPos2) = 1;
end

if ~isempty(pos3)
    result(linPos3) = 2;
end

if ~isempty(pos4)
    result(linPos4) = 3;
end

if ~isempty(pos5)
    result(linPos5) = 4;
end
