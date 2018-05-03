set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set




% dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/clustering/csv/proteome/';
filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/highdose/dmso/*');
filesdmso = filesdmso(3:end);
filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/highdose/soraf/*');
filessoraf = filessoraf(3:end);

dmsodir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/highdose/dmso/';
sorafdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/highdose/soraf/';


for i = 1:size(filesdmso,1)
    dmso = readtable([dmsodir filesdmso(i).name]);
    soraf = readtable([sorafdir filessoraf(i).name]);

    dmso = [dmso(:,10:24) dmso(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
    soraf = [soraf(:,10:24) soraf(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling


    for j = 1:size(dmso,2)
        colname = dmso.Properties.VariableNames(j);
        dmso_col = dmso.(colname{1});
        soraf_col = soraf.(colname{1});
        figure;
        histogram(dmso_col,100,'FaceColor','blue')
        hold on
        histogram(soraf_col,100,'FaceColor','red')
        title(colname{1}, 'Interpreter', 'none')
        xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',14,'fontname','Arial');     
        legend('DMSO','Sorafenib','Location','northwest')
        hold off

        %Save figure 
        fn1 = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/hist/signaling/highdose/';
        fn2 = filesdmso(i).name(28:30);
        fn3 = colname{1}
        fn4 = '_hist.png';
        fn = [fn1 fn2 fn3 fn4];
        print(fn, '-opengl', '-dpng','-r600')

    end
end

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set




% dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/clustering/csv/proteome/';
filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/meddose/dmso/*');
filesdmso = filesdmso(3:end);
filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/meddose/soraf/*');
filessoraf = filessoraf(3:end);

dmsodir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/meddose/dmso/';
sorafdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/meddose/soraf/';


for i = 1:size(filesdmso,1)
    dmso = readtable([dmsodir filesdmso(i).name]);
    soraf = readtable([sorafdir filessoraf(i).name]);

    dmso = [dmso(:,10:24) dmso(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
    soraf = [soraf(:,10:24) soraf(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling


    for j = 1:size(dmso,2)
        colname = dmso.Properties.VariableNames(j);
        dmso_col = dmso.(colname{1});
        soraf_col = soraf.(colname{1});
        figure;
        histogram(dmso_col,100,'FaceColor','blue')
        hold on
        histogram(soraf_col,100,'FaceColor','red')
        title(colname{1}, 'Interpreter', 'none')
        xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',14,'fontname','Arial');     
        legend('DMSO','Sorafenib','Location','northwest')
        hold off

        %Save figure 
        fn1 = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/hist/signaling/meddose/';
        fn2 = filesdmso(i).name(27:29)
        fn3 = colname{1}
        fn4 = '_hist.png';
        fn = [fn1 fn2 fn3 fn4];
        print(fn, '-opengl', '-dpng','-r600')

    end
end

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set




% dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/clustering/csv/proteome/';
filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/lowdose/dmso/*');
filesdmso = filesdmso(3:end);
filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/lowdose/soraf/*');
filessoraf = filessoraf(3:end);

dmsodir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/lowdose/dmso/';
sorafdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/signaling/lowdose/soraf/';


for i = 1:size(filesdmso,1)
    dmso = readtable([dmsodir filesdmso(i).name]);
    soraf = readtable([sorafdir filessoraf(i).name]);

    dmso = [dmso(:,10:24) dmso(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
    soraf = [soraf(:,10:24) soraf(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling


    for j = 1:size(dmso,2)
        colname = dmso.Properties.VariableNames(j);
        dmso_col = dmso.(colname{1});
        soraf_col = soraf.(colname{1});
        figure;
        histogram(dmso_col,100,'FaceColor','blue')
        hold on
        histogram(soraf_col,100,'FaceColor','red')
        title(colname{1}, 'Interpreter', 'none')
        xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',14,'fontname','Arial');     
        legend('DMSO','Sorafenib','Location','northwest')
        hold off

        %Save figure 
        fn1 = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/tfRatios/hist/signaling/lowdose/';
        fn2 = filesdmso(i).name(27:29)
        fn3 = colname{1}
        fn4 = '_hist.png';
        fn = [fn1 fn2 fn3 fn4];
        print(fn, '-opengl', '-dpng','-r600')

    end
end

