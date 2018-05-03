clear;clc;

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set



% % dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/csv/proteome/';
filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/proteome/highdose/dmso/*');
filesdmso = filesdmso(3:end);
filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/proteome/highdose/soraf/*');
filessoraf = filessoraf(3:end);

dmsodir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/proteome/highdose/dmso/';
sorafdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/proteome/highdose/soraf/';

%Here have all soraf for a given dose, timepoints and sig sets are in file
%name. 
%know there will be 4 tps, 2 abs for sig and 3 tps, 3 abs for prot
%split on ab, do timepoints in order, plot histograms in subplot form
% {'highdose_soraf_ProcessedDataP8_proteome3.csv'}

ab_sets = {'proteome1','proteome2','proteome3'};
for i = 1:length(ab_sets)
	matches = strfind({filesdmso(:).name},ab_sets{i});
    filessoraf_absplit{i} = filessoraf(find(~cellfun(@isempty,matches)));
    filesdmso_absplit{i} = filesdmso(find(~cellfun(@isempty,matches)));
end

c = containers.Map;
c('1') = '1 hr';
c('2') = '2 hr';
c('3') = '4 hr';
c('4') = '24 hr';
c('5') = '24 hr';
c('6') = '2 days';
c('7') = '3 days';
c('8') = '5 days';
 
        
        
for ab = 1:length(filessoraf_absplit)   %ab set
    %len, 1=7,2=27,3=28
    if ab < 3
        len = 27;
    else
        len = 28;
    end
    for j = 1:len
         figure('Position', [100, 100, 1750, 750]);
        for i = 1:size(filessoraf_absplit{ab},1) %time point
            dmso = readtable([dmsodir filesdmso_absplit{ab}(i).name]);
            soraf = readtable([sorafdir filessoraf_absplit{ab}(i).name]);
            
            dmso = [dmso(:,9:20) dmso(:,25:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
            soraf = [soraf(:,9:20) soraf(:,25:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling


%         for j = 1:size(dmso,2)  %observable 
            colname = dmso.Properties.VariableNames(j);
            dmso_col = dmso.(colname{1});
            soraf_col = soraf.(colname{1});
            if length(dmso_col) > length(soraf_col)
                soraf_col(end+1:length(dmso_col)) = NaN;
            elseif length(soraf_col) > length(dmso_col)
                dmso_col(end+1:length(soraf_col)) = NaN;
            end
            combined = [dmso_col soraf_col];
            
            
            
            subplot(1,3,i)
            distributionPlot(combined,'color',{'b','r'}); % defaults 
            tit = [colname{1} ' ' c(filesdmso_absplit{ab}(i).name(29))]; %only works for high dose - filename length dependent 
            title(tit, 'Interpreter', 'none')
%             xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
            set(gca,'XTickLabel',{'DMSO','Sorafenib'},'fontsize',18,'fontweight','b','fontname','Arial');
            ylabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
            set(gca,'fontsize',14,'fontname','Arial');     
%             legend('DMSO','Sorafenib','Location','northwest')
            hold off
            

        end
                %Save figure 
        fn1 = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/violinPlots/';
        fn2 = [ab_sets{ab} '/'];
        fn3 = colname{1}
        fn4 = '_violin.png';
        fn = [fn1 fn2 fn3 fn4];
        print(fn, '-opengl', '-dpng','-r600')
    end
end   
    