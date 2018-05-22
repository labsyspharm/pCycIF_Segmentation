function buildViolinPlots()

    warning('off','all')
    set(0,'DefaultFigureVisible','off');
    set(0,'defaultfigurecolor',[1 1 1])
    parpool(10)

    %TO DO:
    %Add signaling dir
    %Automate file/folder handling, make input variable 
    %Change output folder/file name to include plate/ab_set

    drugs = dir('/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/spreadsheets/signaling/highdose/*');

    parfor index = 3:12 %length(drugs) %Cut out ., .., dmso at end

        drugIndex = index;
        drugName = drugs(drugIndex).name;
        filesDrug = dir(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/spreadsheets/signaling/highdose/' drugName '/*']);
        filesDrug = filesDrug(3:end);
        
        %Output drug being processed to console
        printDrug(drugName)
        
        %Keep DMSO separate, need at every iteration. 
        filesdmso = dir('/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/spreadsheets/signaling/highdose/DMSO/*');
        filesdmso = filesdmso(3:end);

        %%%%%%%%%%%%%%%%


        dmsodir = '/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/spreadsheets/signaling/highdose/DMSO/';
        drugdir = ['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/spreadsheets/signaling/highdose/' drugName '/']; 
% 
%         %Here have all soraf for a given dose, timepoints and sig sets are in file
%         %name. 
%         %know there will be 4 tps, 2 abs for sig and 3 tps, 3 abs for prot
%         %split on ab, do timepoints in order, plot histograms in subplot form
% 
        ab_sets = {'signaling1','signaling2'};
        
        filesdrug_absplit = cell(1,2);
        filesdmso_absplit = cell(1,2);

        matches = strfind({filesdmso(:).name},ab_sets{1});
        filesdrug_absplit{1} = filesDrug(find(~cellfun(@isempty,matches)));
        filesdmso_absplit{1} = filesdmso(find(~cellfun(@isempty,matches)));

        matches = strfind({filesdmso(:).name},ab_sets{2});
        filesdrug_absplit{2} = filesDrug(find(~cellfun(@isempty,matches)));
        filesdmso_absplit{2} = filesdmso(find(~cellfun(@isempty,matches)));     
        
        c = containers.Map;
        c('1') = '1 hr';
        c('2') = '2 hr';
        c('3') = '4 hr';
        c('4') = '24 hr';
        c('5') = '24 hr';
        c('6') = '2 days';
        c('7') = '3 days';
        c('8') = '5 days';
            
    for ab = 1:2   %ab sets (2 or 3)
        
        %This block is simply to find the right number of observables for
        %an ab set. There is probabl a more efficient way to do this
        dmso = readtable([dmsodir filesdmso_absplit{ab}(1).name]);
        dmso = dmso(:,3:end); %ignore two indexing columns
        toRemove = strfind(dmso.Properties.VariableNames,'DNA');
        colNamesToRemove = dmso.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
        dmso(:,colNamesToRemove) = [];  %Remove DNA columns
        obsv =  size(dmso.Properties.VariableNames,2);
        %%%%%%%%
        
        for j = 1:obsv   %observables
             figure('Position', [100, 100, 1750, 750]);

%             for fi = 1:size(filesdrug_absplit{ab},1) %time point (files)
            for fi = 1:4 %time point (files)
                dmso = readtable([dmsodir filesdmso_absplit{ab}(fi).name]);
                drug = readtable([drugdir filesdrug_absplit{ab}(fi).name]);  

                %Remove unwanted columns from both dmso and drug tables
                dmso = dmso(:,3:end); %ignore two indexing columns
                toRemove = strfind(dmso.Properties.VariableNames,'DNA');
                colNamesToRemove = dmso.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
                dmso(:,colNamesToRemove) = [];  %Remove DNA columns
                
                drug = drug(:,3:end); %ignore two indexing columns
                toRemove = strfind(drug.Properties.VariableNames,'DNA');
                colNamesToRemove = drug.Properties.VariableNames(find(~cellfun(@isempty,toRemove)));
                drug(:,colNamesToRemove) = [];  %Remove DNA columns
                              
          
                colname = dmso.Properties.VariableNames(j);
                dmso_col = dmso.(colname{1});
                drug_col = drug.(colname{1});
                if length(dmso_col) > length(drug_col)
                    drug_col(end+1:length(dmso_col)) = NaN;
                elseif length(drug_col) > length(dmso_col)
                    dmso_col(end+1:length(drug_col)) = NaN;
                end
                combined = [dmso_col drug_col];

                subplot(1,4,fi)
                distributionPlot(combined,'color',{'b','r'}); % defaults
                %                     tit = [colname{1} ' ' c(filesdmso_absplit{ab}(i).name(67))]; %only works for high dose - filename length dependent
                fn = filesdmso_absplit{ab}(fi).name
                plateCell = strsplit(fn,'_signal')
                plateNum = plateCell{1}(end)
                tit = [colname{1} ' ' c(plateNum)]; %only works for high dose - filename length dependent
                title(tit, 'Interpreter', 'none')
                %             xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
                set(gca,'XTickLabel',{'DMSO',drugName},'fontsize',18,'fontweight','b','fontname','Arial');
                ylabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
                set(gca,'fontsize',14,'fontname','Arial');
                %             legend('DMSO','Sorafenib','Location','northwest')
                
            end
            % %     %
            printDrug(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName])
            if ~exist(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName],'dir')
                mkdir(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName])
            end
            if ~exist(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName '/' ab_sets{ab}],'dir')
                mkdir(['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName '/' ab_sets{ab}])
            end
            
            %Save figure
            fn1 = ['/home/rps21/Cardiomyocite/segmentation/plotResults/tfRatios/allDrugViolinPlots/' drugName '/' ab_sets{ab} '/'];
            fn3 = colname{1}
            fn4 = '_violin.png';
            fn = [fn1 fn3 fn4];
            print(fn, '-dpng','-r300')
        end
    end
    end
    delete(gcp)
end