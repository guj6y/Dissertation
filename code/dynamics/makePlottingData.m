if ~exist('persistences','var')
run('parseData.m')
end

cd ../figures/

numberOfFactors = numel(nFacts); 

 
nFPar = numel(fParAll0);

alpha = 0.05;

speciesCategories = {'all','para','free','basal'};
statisticUsed = {'per','bio','cv','abc-common','abc-natural'};
%make split data ('on' and 'off' for the 2-split; numbered
%for the full split
data = biomasses;

iiHeaderDiffSplit = 'x,y1,y2,y3,y4,m1,m2,m3,m4';
allStats = {persistences,biomasses,cvs,abc.common,abc.natural};
%Think about automating this?
statCode = statisticUsed{1};
stat = persistences;

for statCode = 1:5
    if statCode == 4
        speciesCategories = {'slopes','rSquared'};
    end
    
    stat = allStats{statCode};
    statID = statisticUsed{statCode};
    
    for spCat = speciesCategories
        
        filenamePrefix = sprintf('%s-%s',statID,spCat{:});

        meanDiffs= nan(nFPar,numberOfFactors);
        marginDiffs = nan(nFPar,numberOfFactors);
        
        for ii = 1:numberOfFactors %numberOfFactors will always be 4.
            jjOn = 2;
            jjOff = 1;
            factSize = nFacts(ii);
            selector = cell(numel(size(stat.(spCat{:}))));
            [selector{:}] = deal(':');
            

            selector2 = cell(numel(size(stat.(spCat{:}))));
            [selector2{:}] = deal(':');
            selector2{4} = [jjOff jjOn];
            
            meansSplit    = zeros(nFPar,factSize);
            means2Split    = zeros(nFPar,2);
            marginsSplit  = zeros(nFPar,factSize);
            margins2Split = zeros(nFPar,2);

            iiHeaderYs = 'x,';
            iiHeaderMs = '';
            formatAllSplit = '%.3f,';
            for jj = 1:factSize %originally was 2, but can handle more.
                selector{ii+2} = jj;
                selector2{ii+2} = jj;
                
                datajj  = permute(stat.(spCat{:})(selector{:}),[1,3,4,5,6,2]);
                datajj  = reshape(datajj,[],nFPar);
                
                datajj2 = permute(stat.(spCat{:})(selector2{:}),[1,3,4,5,6,2]);
                datajj2  = reshape(datajj2,[],nFPar);

                if jj == jjOff
                    dataOff = datajj2;
                elseif jj == jjOn
                    dataOn = datajj2;
                end

                meanjj  = mean(datajj,'omitnan');
                stdjj   = std(datajj,'omitnan');
                njj = sum(isfinite(datajj));
                
                meanjj2  = mean(datajj2,'omitnan');
                stdjj2   = std(datajj2,'omitnan');
                njj2 = sum(isfinite(datajj2));

                nObs = sum(isfinite(meanjj));
                nObs2 = sum(isfinite(meanjj2));
                %Need to really think about how to apply Bonferroni here.
                %These are useful if we are doing an all for all bonferroni. 
                m2 = 2*nObs;
                mNFact = factSize*nObs;
                
                tCrit2 = tinv(1-alpha/(2*m2),njj-1);
                tCritN = tinv(1-alpha/(2*mNFact),njj-1);
                
                marginsSplit(:,jj) = tCritN.*stdjj./sqrt(njj);
                margins2Split(:,jj) = tCrit2.*stdjj2./sqrt(njj2);

                meansSplit(:,jj) = meanjj;
                means2Split(:,jj) = meanjj2;

                iiHeaderYs = strcat(iiHeaderYs,sprintf('y%u,',jj));
                iiHeaderMs = strcat(iiHeaderMs,sprintf('m%u,',jj));
                formatAllSplit = strcat(formatAllSplit,'%.9e,%.9e,');
                selector2{4} = [jjOff jjOn];
            end
            dataDiff = dataOn-dataOff;
            nDiff = sum(isfinite(dataDiff));
            
            meanDiff = mean(dataDiff,'omitnan');
            stdDiff = std(dataDiff,'omitnan');
            mDiff = sum(isfinite(meanDiff));
            meanDiffs(:,ii) = meanDiff;
            tCritDiff = tinv(1-alpha/(2*mDiff),nDiff-1);
            marginDiff = stdDiff./sqrt(nDiff).*tCritDiff;
            
            
            marginDiffs(:,ii) = marginDiff;
            
            iiHeaderMs(end) = [];
            formatAllSplit(end) = [];
            
            iiHeaders2Split = 'x,yOn,mOn,yOff,mOff';
            split2FileName = sprintf('%s-2-subplot-%u',filenamePrefix,ii);
            fid2Split = fopen(split2FileName,'w');
%Need to get a format string of %.9f, entries.
            fprintf(fid2Split,'%s',iiHeaders2Split);
            fprintf(fid2Split,'\n%.2f,%.9e,%.9e,%.9e,%.9e',...
                              [fParAll0'...
                              ,means2Split(:,jjOn)...
                              ,margins2Split(:,jjOn)...
                              ,means2Split(:,jjOff)...
                              ,margins2Split(:,jjOff)...
                              ]'...
                   );
            fclose(fid2Split);
            
            iiHeaderFullSplit = strcat(iiHeaderYs,iiHeaderMs);

            splitNFileName = sprintf('%s-full-subplot-%u',filenamePrefix,ii);
            fidNSplit = fopen(splitNFileName,'w');
            fprintf(fidNSplit,'%s',iiHeaderFullSplit);
            fprintf(fidNSplit,sprintf('\n%s',formatAllSplit)...
                              ,[fParAll0'...
                              ,meansSplit...
                              ,marginsSplit]'...
                   );
            fclose(fidNSplit);
            selector{ii+2} = ':';
            selector2{4} = [jjOff jjOn];
        end
        diffFilename = sprintf('%s-diffs',filenamePrefix);
        fidDiff = fopen(diffFilename,'w');
        diffFormat = strcat('%.3f,',repmat('%.9e,',1,7),'%.9e');
        fprintf(fidDiff,'%s',iiHeaderDiffSplit);
        fprintf(fidDiff,sprintf('\n%s',diffFormat),...
                        [fParAll0'...
                        ,meanDiffs...
                        ,marginDiffs...
                        ]'...
                );
        fclose(fidDiff);
    end
end
