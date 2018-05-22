function [] = histogramsMarkerVsControl(controlCol)

set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:


files = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output/NCRatio/*.txt');   %Make base an input variable
files = files(3:end);
filesdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output/NCRatio/';

%Gather control files and treatment files, based on control being a given
%column. Eventually change to be an input variable. 
if ~controlCol
    controlCol = '11';
end

fileFields = fieldnames(files)';
fileFields{2,1} = {};
controlFiles = struct(fileFields{:});   %Should try to preallocate these if possible
treatmentFiles = struct(fileFields{:});

for i = 1:size(files,1) %time point
    fileName = files(i).name;
    if contains(fileName,controlCol)
        controlFiles=[controlFiles,files(i)];
    else
        treatmentFiles=[treatmentFiles,files(i)];
    end  
end
    

%Loop through control
%loop through treatment, checking for matching row
    %trick will be extracting row info
        %split = strsplit(files(i).name,'_')
        %split{2}(1) 
        %gives letter, but feels fragile
%figure - plot both together
for i = 1:size(controlFiles,2) 
    controlData = readtable([filesdir controlFiles(i).name]);
    colNames = controlData.Properties.VariableNames;
    colNamesToUse = colNames(~startsWith(colNames,{'Bleach','Var','Unnamed'}));
    controlData = controlData(:,colNamesToUse) ;
    colRow = strsplit(controlFiles(i).name,'_');
    rowToUse = colRow{2}(1);     %gives letter, but feels fragile
    
    for j = 1:size(treatmentFiles,2)
        colRow = strsplit(treatmentFiles(j).name,'_');
        row = colRow{2}(1); %gives letter, but feels fragile
        
        if strcmp(row,rowToUse)
           
            treatmentData = readtable([filesdir treatmentFiles(j).name]);
            treatmentData = treatmentData(:,colNamesToUse) ;
    
            for k = 1:length(colNamesToUse)  %observable
                colname = colNamesToUse{k};
%                 colname = fileData.Properties.VariableNames(j);
                markerControl = controlData.(colname);
                markerTreatment = treatmentData.(colname);              
                                
                figure('Position', [100, 100, 1000, 750]);        
                histogram(markerControl,'FaceColor','blue')
                hold on
                histogram(markerTreatment,'FaceColor','red')
                title(colname, 'Interpreter', 'none')
                xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
                set(gca,'fontsize',14,'fontname','Arial');
                legend('Control','Treatment','Location','northwest')
                hold off

                folderBase = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/processingAndVisualization/testOutput/';  %Make input variable or extrapolated from above
                inputFile = strsplit(treatmentFiles(i).name,'.');
                inputFileName = [inputFile{1} '/'];
                %check/create folder here
                if ~exist([folderBase inputFileName],'dir')
                    mkdir([folderBase inputFileName])
                end
                marker = colname %Create better output for tracking
                fileEnd = '_hist.png';
                fn = [folderBase inputFileName marker fileEnd];
                print(fn, '-opengl', '-dpng','-r300')

            end
        end
    end
end
