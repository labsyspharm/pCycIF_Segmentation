set(0,'DefaultFigureVisible','off');
set(0,'defaultfigurecolor',[1 1 1])
%TO DO:
%Fix file name inputs

files = dir('/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output/NCRatio/*.txt');   %Make base an input variable
files = files(3:end);

%Can we extract row/col/fld from here?
%loop through and find files for control when specified here?

filesdir = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/output/NCRatio/';
  
for i = 1:size(files,1) %time point
    fileData = readtable([filesdir files(i).name]);
    colNames = fileData.Properties.VariableNames;
    colNamesToUse = colNames(~startsWith(colNames,{'Bleach','Var','Unnamed'}));
    fileData = fileData(:,colNamesToUse) ;
    
    for j = 1:size(colNamesToUse,2)  %observable
        colname = fileData.Properties.VariableNames(j);
        marker_col = fileData.(colname{1});
        
        figure('Position', [100, 100, 1000, 750]);        
        histogram(marker_col,'FaceColor','blue')
        title(colname{1}, 'Interpreter', 'none')
        xlabel('log2 Intensity','fontsize',18,'fontweight','b','fontname','Arial');
        set(gca,'fontsize',14,'fontname','Arial');
        hold off
    
        folderBase = '/home/bobby/Dropbox/MATLAB/cardiotoxCycif/segmentation/pCycIF_Segmentation/processingAndVisualization/testOutput/allMarkers/';  %Make input variable or extrapolated from above
        inputFile = strsplit(files(i).name,'.');
        inputFileName = [inputFile{1} '/'];
        %check/create folder here
        if ~exist([folderBase inputFileName],'dir')
            mkdir([folderBase inputFileName])
        end
        marker = colname{1} %Create better output for tracking
        fileEnd = '_hist.png';
        fn = [folderBase inputFileName marker fileEnd];
        print(fn, '-opengl', '-dpng','-r300')

    end
    
end
