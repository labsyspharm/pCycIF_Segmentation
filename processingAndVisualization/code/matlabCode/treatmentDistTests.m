%looping through folders, files?
%sorafProcessedDataP1signaling1 2011x43
ctrl = dmsoProcessedDataP1signaling1;
treatment = sorafProcessedDataP1signaling1;
for i = 1:size(ctrl,2)
    i
    [h,p] = ttest2(ctrl(:,i),treatment(:,i));
    [h2,p2] = kstest2(ctrl(:,i),treatment(:,i));
    [p3,h3] = ranksum(ctrl(:,i),treatment(:,i));
    hTot(1,i) = h;
    hTot(2,i) = h2;
    hTot(3,i) = h3;
    pTot(1,i) = p;
    pTot(2,i) = p2;
    pTot(3,i) = p3;
end