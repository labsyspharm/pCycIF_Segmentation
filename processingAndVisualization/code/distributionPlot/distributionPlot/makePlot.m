% r = rand(1000,1); 
% rn = randn(1000,1)*0.38+0.5; 
% rn2 = [randn(500,1)*0.1+0.27;randn(500,1)*0.1+0.73]; 
%     rn2=min(rn2,1);rn2=max(rn2,0); 
% 
% 
% figure 
% ah(1)=subplot(2,4,1:2); 
% boxplot([r,rn,rn2]) 
% 
% figure;
% ah(2)=subplot(2,4,3:4); 
% distributionPlot([r,rn,rn2],'histOpt',2); % histOpt=2 works better for uniform distributions than the default 
% set(ah,'ylim',[-1 2]) 
% 
% %--additional options 
% figure;
% data = [randn(100,1);randn(50,1)+4;randn(25,1)+8]; 
% subplot(2,4,5) 
% distributionPlot(data); % defaults 
% 
% figure;
% subplot(2,4,6) 
% distributionPlot(data,'colormap',copper,'showMM',5,'variableWidth',false) % show density via custom colormap only, show mean/std, 
% 
% figure;
% subplot(2,4,7:8) 
% distributionPlot({data(1:5:end),repmat(data,2,1)},'addSpread',true,'showMM',false,'histOpt',2) %auto-binwidth depends on # of datapoints; for small n, plotting the data is useful
% 


%--additional options 
figure;
subplot(1,8,1) 
distributionPlot(dmso_col,'color','b'); % defaults 
subplot(1,8,2) 
axis off
distributionPlot(soraf_col,'color','r'); % defaults 


% dmsocoltest(end+1:847) = NaN
%find which is longer
%add end+1:length() NaN
figure;
combined = [dmsocoltest soraf_col];
distributionPlot(combined,'color',{'b','r'}); % defaults 




%       color : uniform coloring of histograms. Supply either a color
%           string ('r'), or a truecolor vector ([1 0 0]). Use a
%           cell array of length nData to specify one color per
%           distribution. Default: 'k'