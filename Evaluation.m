function [P,R,F1] = Evaluation(TP,FP,FN)

if TP==0 & FP==0
    P = 0;
else
    P = TP/(TP+FP);
end

if TP==0 & FN==0
    P = 0;
else
    R = TP/(TP+FN);
end

if P==0 & R==0
    F1 = 0;
else
    F1 = 2*((P*R)/(P+R));
end

end