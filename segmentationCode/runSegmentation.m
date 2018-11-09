function [] = runSegmentation_Caitlin(parametersFile)
% function CaitlinpCycifNucCytoSegmentationCPU(parentPath,modelPath,modelCatPath,FFCPath,varargin)
%this function requires input of nuclei stack range. It assumes that every
%stack beyond that to the end is a cytoplasmic stain. Marker controlled
%watershed based on distance transform of nuclei channel is employed to
%separate nuclei clumps.

javaaddpath('matlabDependencies/bfmatlab/bioformats_package.jar')
javaaddpath('matlabDependencies/yaml/java/snakeyaml-1.9.jar')

[params] = YAML.read(parametersFile);
NucMaskChan = str2num(params.NucMaskChan{1});
CytoMaskChan = str2num(params.CytoMaskChan{1});

parentPath = params.parentPath;
analysisPath = params.outputPath;
modelPath = params.modelPath;   %.mat, static
modelCatPath = params.modelCatPath; %.mat, static
FFCPath = params.FFCPath;   %Tif, static
% javaPath = params.javaPath;
numCycles = params.numCycles;
Row = params.row;
Col = params.col;
saveFig = params.saveFig;
cytoMethod = params.cytoMethod;
MedianIntensity = params.MedianIntensity;
saveMasks = params.saveMasks;
applyFFC = params.applyFFC;
useRFNuc = params.useRFNuc;
segmentCytoplasm = params.segmentCytoplasm;

%Set NucMaskChan from numCycles = 5;
NucMaskChan = 2:numCycles;

%% initialization
drive='Y';
% parentPath = [drive ':\sorger\data\IN Cell Analyzer 6000\Connor\Fixed MCF10 Common\20x full exp\20180905_Updated'];
folders=dir([parentPath filesep '*Plate*']);

% FFCPath = [drive ':\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\Real experiment\All Cyles\FFC'];
switch applyFFC
    case 'both'
        dfp = volumeRead([FFCPath filesep 'allDF.tif']);
        ffp = volumeRead([FFCPath filesep 'allFFP.tif']);
    case 'ffonly'
        ffp = volumeRead([FFCPath filesep 'allFFPOnly.tif']);
    case 'none'
end

% modelPath = [drive ':\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\nucleiTrainingSet\halfsize'];
%pragma necessary for Matlab compiler to recognize treeBagger class when
%loading from .mat file
%#function treeBagger
load ([modelPath filesep 'model.mat'])
% modelCatPath = [drive ':\sorger\data\IN Cell Analyzer 6000\Connor\MCF10 Common\20x full exp\cateninTrainingSet'];
load ([modelCatPath filesep 'modelCatManual.mat'])


%%
for iFolder = 1:numel(folders)
    testPath = [parentPath filesep folders(iFolder).name];
   
    if exist([testPath filesep 'analysis'],'dir')~=7
        mkdir([testPath filesep 'analysis'])
        
    end
    %     analysisPath = [testPath filesep 'analysis'];
    
    row = Row(1):Row(2);
    col = Col(1):Col(2);
    for iRow = 1:numel(row)
        for iCol = 1:numel(col)
            for iField = 1:9
                files = dir([testPath filesep  row(iRow) sprintf('%.2d', col(iCol)) '_fld' int2str(iField) '*.tif']);
                for iFile = 1:numel(files)
               
                    tic
                    fileName = files(iFile).name;
                    [~,name,ext] = fileparts(fileName) ;
                    
                    I = volumeRead([testPath filesep fileName]);
                    nucleiStack = [1 size(I,3)/4];
                    nucleiMaskChan = NucMaskChan;
                    cytoChanRange = (size(I,3)-nucleiStack(2));
                    cytoChanStart =size(I,3)/4+1;
                    cytoChanEnd = size(I,3);
                    
                    nucleiImage = I(:,:,nucleiMaskChan(1):nucleiMaskChan(2));
                    nucleiImage = max(nucleiImage,[],3);
                    nucleiImage = double(nucleiImage);
                    nucleiImage = nucleiImage/65535;%max(nucleiImage(:));
                    
                    if (useRFNuc)
                        %% apply random forest model to generate classProbs
                        F = pcImageFeatures(imresize(nucleiImage,0.5,'nearest'),model.sigmas,model.offsets,model.osSigma,model.radii,...
                            model.cfSigma,model.logSigmas,model.sfSigmas,model.ridgeSigmas,model.ridgenangs,...
                            model.edgeSigmas,model.edgenangs,model.nhoodEntropy,model.nhoodStd);
                        [imL,classProbs] = imClassify(F,model.treeBag,100);
                        nucleiClass=classProbs(:,:,3);
                        
                        %% markers based on log filter on classProbs 3
                        logNuclei=  imgaussfilt3(filterLoG(nucleiClass,4),2.5); %3.5, 2
                        logfgm = imregionalmax(logNuclei);
                        logfgm = (logNuclei>0).*(logfgm==1);
                        %                           figure,imshowpair(logfgm,sqrt(normalize(imresize(nucleiClass,1))))
                        
                        
                        %% contours for watershed
                        nucleiBlur = medfilt2(nucleiImage,[11 11]);
                        hy = fspecial('sobel');
                        hx = hy';
                        Iy = imfilter(double(nucleiBlur), hy, 'replicate');
                        Ix = imfilter(double(nucleiBlur), hx, 'replicate');
                        gradmag = sqrt(Ix.^2 + Iy.^2);
                        contours = imresize(classProbs(:,:,2),[size(gradmag,1) size(gradmag,2)]).*gradmag;
                        %% apply watershed transform
                        bgm = imerode(classProbs(:,:,1)>0.3,strel('disk',2));
                        gradmag2= imimposemin(contours,imresize(bgm|logfgm,[size(contours,1) size(contours,2)]));
                        foregroundMask= watershed(gradmag2);
                        
                        %% process mask
                        singleNuclei = bwareafilt(foregroundMask>0,[0 799]);
                        multiNuclei  = imresize(imfill(bwareafilt(foregroundMask>0,[800 1000]),'holes'),2);
                        
                        
                        %% marker controlled watershed
                        if sum(sum(multiNuclei))>0
                            imgDist =-bwdist(~multiNuclei);
                            imgauss=imhmin(imgDist,1);
                            Imax = imregionalmin(imgauss);
                            imgDist=imimposemin(imresize(contours,[size(Imax,1) size(Imax,2)]),Imax);
                            imgDist(~multiNuclei)=-inf;
                            imgLabel=watershed(imgDist);
                            imgLabel = bwareafilt(imgLabel>0,[0 5000]);
                        else
                            imgLabel =multiNuclei;
                        end
                        
                        
                        allNuclei = imresize(singleNuclei,2)| imgLabel;
                        allNuclei = bwareaopen(allNuclei,500); %500
                        stats=regionprops(bwlabel(allNuclei),imresize(classProbs(:,:,3),[size(allNuclei,1) size(allNuclei,2)]),'MeanIntensity','Area');
                        idx = find([stats.MeanIntensity] > 0.2);
                        nucleiMask = ismember(bwlabel(allNuclei),idx);
                        
                        statsNM=regionprops(imresize(nucleiMask,[size(allNuclei,1) size(allNuclei,2)]),'Area');
                        largestNucleusArea=prctile([statsNM.Area],95);
                        
                        %                        figure,imshowpair(edge(nucleiMask),imresize(nucleiImage,2))
                        
                    else
                        %%use presaved label masks
                        nucleiMask = imread([testPath filesep name '_nucleiLM' ext]);
                        %                         %%use probability maps from UNet
                        %                         nucProbMap = volumeRead([testPath filesep name '_PM.tif']);
                        %                         radii = [4 6];
                        %                         count=0;
                        %                         A=[];
                        %                         for iradii = radii
                        %                         count=count+1;
                        %                         [A(:,:,count),K] = circcentliklK(nucProbMap(:,:,2),iradii,1.5,8,true);
                        %                         end
                        %
                        %                         fgm = max(A,[],3)>400;imshowpair(fgm,nucProbMap(:,:,2))
                        %
                        %                         gradmag = imimposemin(nucProbMap(:,:,2),fgm);
                        %                         nuclei = watershed(gradmag);
                        % %                         imshowpair(nuclei==0,sqrt(normalize(nucProbMap(:,:,1))))
                        %
                        %                         label =bwconncomp(nuclei>0);
                        %                         stats=regionprops(label,imtophat(nucProbMap(:,:,1),strel('disk',15)),'MeanIntensity');
                        %                         idx=find([stats.MeanIntensity]>thresholdOtsu(imtophat(nucProbMap(:,:,1),strel('disk',10))));
                        %                         nucleiMask = ismember(labelmatrix(label), idx) > 0;
                        %                         nucleiMask = imresize(bwareaopen(nucleiMask,10),4);
                    end
                    
                    %% cytoplasm segmentation
                    switch segmentCytoplasm
                        case 'segmentCytoplasm'
                            cyto = max(I(:,:,CytoMaskChan(1):CytoMaskChan(2)),[],3);
                            switch cytoMethod
                                case 'RF'
                                    F = pcImageFeatures(imresize(double(cyto)/65335,0.5,'bilinear'),modelCat.sigmas,modelCat.offsets,...
                                        modelCat.osSigma,modelCat.radii,modelCat.cfSigma,modelCat.logSigmas,modelCat.sfSigmas,...
                                        modelCat.ridgeSigmas,modelCat.ridgenangs,modelCat.edgeSigmas,modelCat.edgenangs,modelCat.nhoodEntropy,...
                                        modelCat.nhoodStd);
                                    [imL,catClassProbs] = imClassify(F,modelCat.treeBag,100);
                                    contours = imresize(catClassProbs(:,:,2),2);
                                    %                             bgm = classProbs(:,:,1)>0.9;
                                    %                             bgm = imresize(bwmorph(bgm,'shrink',Inf),4);
                                    bgm =imresize(bwmorph( imgaussfilt3(cyto,2)<100,'thin',Inf),2);
                                    cytograd= imimposemin(imresize(contours,[size(nucleiMask,1) size(nucleiMask,2)]),bgm|nucleiMask);
                                    cellMask= watershed(cytograd);
                                    
                                case 'distanceTransform'
                                    bgLoG=imgaussfilt3(filterLoG(max(cyto(:))-cyto,10),10);
                                    bgMax = imregionalmax(bgLoG);
                                    threshold = median(bgLoG(bgMax==1));
                                    bgm=imresize(uint16(bgMax.*(bgLoG<threshold)),2);
                                    contours =-bwdist(imresize(cyto<thresholdOtsu(cyto),2));
                                    cytograd= imimposemin(imresize(contours,[size(nucleiMask,1) size(nucleiMask,2)]),bgm|nucleiMask);
                                    cellMask=watershed(cytograd);
                                    thresholdMask=mean(cyto(bgMax==1));
                                    preCellMask = bwareaopen(imgaussfilt3(cyto,1)>thresholdMask,400);
                                    cellMask = cellMask.*cast(imresize(preCellMask,2),class(cellMask));
                                    cellMask = bwlabel(bwareaopen(cellMask,round(largestNucleusArea)));
                                    
                                case 'contours'
                                    
                                    contours = normalize(steerableDetector(im2double(cyto),2,1.5));
                                    
                                case 'ring'
                                    
                            end
                            
                            nucleiMask = cast(nucleiMask,class(cellMask)).*cellMask;
                            
                            %clean this up!
                            test=cast(~ismember(unique(cellMask),unique(nucleiMask)),class(unique(cellMask))).*unique(cellMask);
                            bgCells=find(test>0);
                            for iTest=bgCells'
                                cellMask(cellMask ==test(iTest))=0;
                            end
                            
                            
                            %% eliminate border cells
                            inCells = imclearborder(cellMask>0);
                            borderCells = cellMask.*cast((cellMask>0)-inCells,class(cellMask));
                            borderIdx = unique(borderCells);
                            nucleiMask_border = ~ismember(nucleiMask,borderIdx);
                            
                            cellMask = bwlabel(inCells);
                            nucleiMask = nucleiMask_border.*cellMask;
                            cytoplasmMask = cellMask - nucleiMask;
                            %                 imshowLinkedTuple(imresize(normalize(im2double(cyto(:,:,1))),2) + im2double(nucleiMask>0),cellMask==0)
                            
                        case 'loadMask'
                            cytoplasmMask = imread([testPath filesep name '_cytoLM' ext]);
                        case 'ignoreCytoplasm'
                            cytoplasmMask = nucleiMask.*0;
                    end
                    %                         imshowLinkedTuple(cellMask==0,imresize(sqrt(normalize(cyto)),2))
                    %% apply flatfield correction using flatfield & darkfield, just flatfield, or no correction
                    switch applyFFC
                        
                        case 'both'
                            for iffp = 1:4
                                for i =1:(size(I,3)/size(ffp,3))
                                    FFCI(:,:,(iffp-1)*size(I,3)/size(ffp,3)+i) = (I(:,:,(iffp-1)*size(I,3)/size(ffp,3)+i)-dfp(:,:,iffp))./ffp(:,:,iffp) ...
                                        *mean(reshape(ffp(:,:,iffp),[1 size(ffp,1)*size(ffp,2)]));
                                end
                            end
                            cytoResized=imresize(FFCI,[size(nucleiMask,1) size(nucleiMask,2)]);
                            
                        case 'ffonly'
                            for iffp = 1:4
                                for i =1:(size(I,3)/size(ffp,3))
                                    FFCI(:,:,(iffp-1)*size(I,3)/size(ffp,3)+i) = (I(:,:,(iffp-1)*size(I,3)/size(ffp,3)+i))./ffp(:,:,iffp) ...
                                        *mean(reshape(ffp(:,:,iffp),[1 size(ffp,1)*size(ffp,2)]));
                                end
                            end
                            cytoResized=imresize(FFCI,[size(nucleiMask,1) size(nucleiMask,2)]);
                            
                        case 'none'
                            cytoResized= imresize(I,[size(nucleiMask,1) size(nucleiMask,2)]);
                    end
                    
                    %% measure intensities from regions
                    meanIntNucTable = zeros(max(nucleiMask(:)),size(I,3));
                    meanIntCytoTable = zeros(max(nucleiMask(:)),size(I,3));
                    medianIntNucTable = zeros(max(nucleiMask(:)),size(I,3));
                    medianIntCytoTable = zeros(max(nucleiMask(:)),size(I,3));
                    centroidCellTable = zeros(max(nucleiMask(:)),2);
                    
                    for iChan = 1: cytoChanEnd
                        
                        nucleiStats=regionprops(nucleiMask,cytoResized(:,:,iChan),'MeanIntensity','Centroid','Area','PixelIdxList');
                        cytoStats=regionprops(cytoplasmMask,cytoResized(:,:,iChan),'MeanIntensity','Centroid','Area','PixelIdxList');
                        
                        if MedianIntensity ==1
                            cytoSlice=cytoResized(:,:,iChan);
                            for iCell = 1: numel(nucleiStats)
                                medianIntNucTable(iCell,iChan) = median(cytoSlice(nucleiStats(iCell).PixelIdxList));
                                medianIntCytoTable(iCell,iChan) = median(cytoSlice(cytoStats(iCell).PixelIdxList));
                            end
                        end
                        
                        meanIntNucTable(:,iChan) = [nucleiStats.MeanIntensity]';
                        meanIntCytoTable(:,iChan) = [cytoStats.MeanIntensity]';
                        
                    end
                    
                    meanIntTable = [meanIntNucTable meanIntCytoTable];
                    medianIntTable = [medianIntNucTable medianIntCytoTable];
                    areaTable = [cat(1,nucleiStats.Area) cat(1,cytoStats.Area)  ];
                    centroidCellTable = cat(1,nucleiStats.Centroid);
                    
                    %% write results to txt file
                    if ~isempty(areaTable)
                        variableNucNamesMeanIntensity = {};
                        variableCytoNamesMeanIntensity = {};
                        variableNucNamesMedianIntensity = {};
                        variableCytoNamesMedianIntensity = {};
                        
                        for ivarName = 1:size(meanIntTable,2)/2
                            variableNucNamesMeanIntensity = cat(2,variableNucNamesMeanIntensity,{['NucleiChannelMeanIntensity' int2str(ivarName)]});
                            variableCytoNamesMeanIntensity = cat(2,variableCytoNamesMeanIntensity,{['CytoplasmChannelMeanIntensity' int2str(ivarName)]});
                            variableNucNamesMedianIntensity = cat(2,variableNucNamesMedianIntensity,{['NucleiChannelMedianIntensity' int2str(ivarName)]});
                            variableCytoNamesMedianIntensity = cat(2,variableCytoNamesMedianIntensity,{['CytoplasmChannelMedianIntensity' int2str(ivarName)]});
                        end
                        writetable(array2table([meanIntTable medianIntTable areaTable centroidCellTable],'VariableNames',[variableNucNamesMeanIntensity variableCytoNamesMeanIntensity...
                            variableNucNamesMedianIntensity variableCytoNamesMedianIntensity 'NucleusArea' 'CytoplasmArea' 'CellPosition_X' 'CellPosition_Y']),...
                            [analysisPath filesep '_' name '_cytoMasked.txt'],'Delimiter','\t')
                    end
                    
                    %% display
                    
                    if saveFig==true
                        
                        % save mask overlay to .fig
                        figure,
                        axs=[];
                        axs = [axs subplot(1,2,1)];  imshow(sqrt(cytoResized(:,:,nucleiMaskChan(1))/max(max(cytoResized(:,:,nucleiMaskChan(1))))))
                        hold on
                        visboundaries(bwboundaries(nucleiMask),'LineWidth',1)
                        
                        for i =1:max(nucleiMask(:))
                            text (nucleiStats(i).Centroid(1),nucleiStats(i).Centroid(2),int2str(i),'Color' ,'r')
                        end
                        
                        axs = [axs subplot(1,2,2)];imshow(sqrt(max(cytoResized(:,:,CytoMaskChan(1):CytoMaskChan(2)),[],3)/max(max(max(cytoResized(:,:,CytoMaskChan(1):CytoMaskChan(2)))))))
                        hold on
                        visboundaries(bwboundaries(cytoplasmMask),'LineWidth',1)
                        
                        for i =1:max(cytoplasmMask(:))
                            text (cytoStats(i).Centroid(1),cytoStats(i).Centroid(2),int2str(i),'Color' ,'r')
                        end
                        linkaxes(axs,'xy')
                        savefig ([analysisPath filesep name '_cytoMasked.fig' ])
                        
                        %save raw upsampled image
                        tiffwriteimj(uint16(imresize(I,[size(nucleiMask,1) size(nucleiMask,2)])), [analysisPath filesep name '_RAW' ext])
                        
                        %save edge mask on raw image
                        unmaskedImage = I/65535;
                        allChan=[];
                        allEdge = imresize(imdilate(edge(cytoplasmMask>0,'Sobel'),strel('square',2)),[size(I,1) size(I,2)]);
                        for i = 1:size(unmaskedImage,3)
                            allChan(:,:,i)  = 255*(im2double(allEdge))+255*(unmaskedImage(:,:,i));
                        end
                        tiffwriteimj(uint8(allChan), [analysisPath filesep name '_Masked' ext])
                        
                        %save FFC upsampled image
                        if  strcmp(applyFFC,'ffonly') ||  strcmp(applyFFC,'both')
                            tiffwriteimj(uint16(imresize(cytoResized,[size(nucleiMask,1) size(nucleiMask,2)])), [analysisPath filesep name '_FFC' ext])
                        end
                        
                        %save nuclear and cytoplasm masks
                        if saveMasks
                            tiffwriteimj(uint16(nucleiMask), [analysisPath filesep name '_nucleiLM' ext])
                            tiffwriteimj(uint16(cytoplasmMask), [analysisPath filesep name '_cytoLM' ext])
                        end
                        
                        close all
                        disp(['Completed ' fileName])
                    end
                    toc
                end
            end
            disp(['Completed well ' row(iRow) num2str(col(iCol))])
            
        end
    end
    
end

end