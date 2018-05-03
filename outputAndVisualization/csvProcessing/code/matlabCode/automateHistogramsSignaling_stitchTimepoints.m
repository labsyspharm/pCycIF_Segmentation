clear;clc;

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set




% % dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/csv/proteome/';
filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/signaling/highdose/dmso/*');
filesdmso = filesdmso(3:end);
filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/signaling/highdose/soraf/*');
filessoraf = filessoraf(3:end);

dmsodir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/signaling/highdose/dmso/';
sorafdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/spreadsheets/signaling/highdose/soraf/';

%Here have all soraf for a given dose, timepoints and sig sets are in file
%name. 
%know there will be 4 tps, 2 abs for sig and 3 tps, 3 abs for prot
%split on ab, do timepoints in order, plot histograms in subplot form
% {'highdose_soraf_ProcessedDataP8_proteome3.csv'}

ab_sets = {'signaling1','signaling2'};
for i = 1:length(ab_sets)
	matches = strfind({filesdmso(:).name},ab_sets{i});
    filessoraf_absplit{i} = filessoraf(find(~cellfun(@isempty,matches)));
    filesdmso_absplit{i} = filesdmso(find(~cellfun(@isempty,matches)));
end

   
        
        
for ab = 1%1:length(filessoraf_absplit)   %ab set
    
    for j = 1:37
         figure('Position', [100, 100, 1750, 750]);
        for i = 1:size(filessoraf_absplit{ab},1) %time point
            dmso = readtable([dmsodir filesdmso_absplit{ab}(i).name]);
            soraf = readtable([sorafdir filessoraf_absplit{ab}(i).name]);

            dmso = [dmso(:,10:24) dmso(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
            soraf = [soraf(:,10:24) soraf(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling


%         for j = 1:size(dmso,2)  %observable 
            colname = dmso.Properties.VariableNames(j);
            dmso_col = dmso.(colname{1});
            soraf_col = soraf.(colname{1});
           
            subplot(1,4,i)
            histogram(dmso_col,50,'FaceColor','blue')
            hold on
            histogram(soraf_col,50,'FaceColor','red')
            title(colname{1}, 'Interpreter', 'none')
            xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
            set(gca,'fontsize',14,'fontname','Arial');     
            legend('DMSO','Sorafenib','Location','northwest')
            hold off

        end
                %Save figure 
        fn1 = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/timecourseHistograms/';
        fn2 = [ab_sets{ab} '/'];
        fn3 = colname{1}
        fn4 = '_hist.png';
        fn = [fn1 fn2 fn3 fn4];
        print(fn, '-opengl', '-dpng','-r600')
    end
end   
    