clear;clc;

set(0,'DefaultFigureVisible','on');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Add signaling dir
%Change output folder/file name to include plate/ab_set

% conditions = {'pRb_high','pRb_low'};
conditions = {'pRb_low'};

for thing = 1:length(conditions)
    drugs = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/*');
    for i = 11:length(drugs) %Cut out ., .., dmso at end
        drugName = drugs(i).name;
        filesDrug = dir(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/' drugName '/' conditions{thing} '/*']);
        filesDrug = filesDrug(3:end);


        %Keep DMSO separate, need at every iteration. 
        filesdmso = dir(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/DMSO/' conditions{thing} '/*']);
        filesdmso = filesdmso(3:end);

        %%%%%%%%%%%%%%%%
        % % dirname = 'home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/csv/proteome/';
        % filesdmso = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/spreadsheets/signaling/highdose/dmso/*');
        % filesdmso = filesdmso(3:end);
        % filessoraf = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/spreadsheets/signaling/highdose/soraf/*');
        % filessoraf = filessoraf(3:end);

        dmsodir = ['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/DMSO/' conditions{thing} '/'];
        drugdir = ['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/pRbSplit/spreadsheets/highdose/' drugName '/' conditions{thing} '/']; 

        %Here have all soraf for a given dose, timepoints and sig sets are in file
        %name. 
        %know there will be 4 tps, 2 abs for sig and 3 tps, 3 abs for prot
        %split on ab, do timepoints in order, plot histograms in subplot form
        % {'highdose_soraf_ProcessedDataP8_proteome3.csv'}

        ab_sets = {'signaling1','signaling2','proteome1','proteome2','proteome3'};
        for i = 1:length(ab_sets)
            matches = strfind({filesdmso(:).name},ab_sets{i});
            filesdrug_absplit{i} = filesDrug(find(~cellfun(@isempty,matches)));
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



        for ab = 1%:length(filesdrug_absplit)   %ab sets

            if ~isempty(filesdmso_absplit{ab})

                %This block is simply to find the right number of observables for
                %an ab set. There is probabl a more efficient way to do this
                dmso = readtable([dmsodir filesdmso_absplit{ab}(1).name]);
                dmso = dmso(:,3:end); %ignore two indexing columns
                toRemove = strfind(dmso.Properties.VariableNames,'DNA');
                colNamesToRemove = dmso.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
                dmso(:,colNamesToRemove) = [];  %Remove DNA columns
                obsv =  size(dmso.Properties.VariableNames,2);
                %%%%%%%%

                for j = 1%:obsv   %observables
                     figure('Position', [100, 100, 1750, 750]);
                    for i = 1:size(filesdrug_absplit{ab},1) %time point
                        dmso = readtable([dmsodir filesdmso_absplit{ab}(i).name]);
                        drug = readtable([drugdir filesdrug_absplit{ab}(i).name]);

        %                 dmso = [dmso(:,10:24) dmso(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
        %                 drug = [drug(:,10:24) drug(:,30:end)]; %Cut indices and DNA columns - needs to be different for proteome and signaling
                        dmso = dmso(:,3:end); %ignore two indexing columns
                        toRemove = strfind(dmso.Properties.VariableNames,'DNA');
                        colNamesToRemove = dmso.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
                        dmso(:,colNamesToRemove) = [];  %Remove DNA columns

                        drug = drug(:,3:end); %ignore two indexing columns
                        toRemove = strfind(drug.Properties.VariableNames,'DNA');
                        colNamesToRemove = drug.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
                        drug(:,colNamesToRemove) = [];  %Remove DNA columns

            %         for j = 1:size(dmso,2)  %observable 
                        colname = dmso.Properties.VariableNames(j);
                        dmso_col = dmso.(colname{1});
                        drug_col = drug.(colname{1});
                        if length(dmso_col) > length(drug_col)
                            drug_col(end+1:length(dmso_col)) = NaN;
                        elseif length(drug_col) > length(dmso_col)
                            dmso_col(end+1:length(drug_col)) = NaN;
                        end
                        combined = [dmso_col drug_col];



%                         subplot(1,3,i)
%                         distributionPlot(combined,'color',{'b','r'}); % defaults 
%                         tit = [colname{1} ' ' c(filesdmso_absplit{ab}(i).name(end-4))]; %only works for high dose - filename length dependent 
%                         title(tit, 'Interpreter', 'none')
%             %             xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
%                         set(gca,'XTickLabel',{'DMSO',drugName},'fontsize',18,'fontweight','b','fontname','Arial');
%                         ylabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
%                         set(gca,'fontsize',14,'fontname','Arial');     
%             %             legend('DMSO','Sorafenib','Location','northwest')
%                         hold off


                    end

%                     if ~exist(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/violinPlots/' drugName],'dir')
%                         mkdir(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/violinPlots/' drugName])
%                     end
%                     if ~exist(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/violinPlots/' drugName '/' conditions{thing} '/'],'dir')
%                         mkdir(['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/violinPlots/' drugName '/' conditions{thing} '/'])
%                     end
% 
%                     %Save figure 
%                     fn1 = ['/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/newSegmentationData/dapiAdded/csvProcessing/nuc_pRbSplit/violinPlots/' drugName '/' conditions{thing} '/'];
%                     fn2 = colname{1}
%                     fn3 = ['_' ab_sets{ab}]
%                     fn4 = '_violin.png';
%                     fn = [fn1 fn2 fn3 fn4];
%                     print(fn, '-opengl', '-dpng','-r600')
                end
        %         ab


            end
        end   
    end
end
    