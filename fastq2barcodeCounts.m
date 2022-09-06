%% this script identify ASKA barcode within reads in multiple fastq files. 
% the script can use two alternative methods:
% 1. Identify an EXACT mantch of a 20bp barcode that immedietely after a 
% PS1 anchor (a conserved area upstream of the barcode). minute per screen 
% 2. Serach for matches to any barcode (15-25bp) in the approx region of
% barcode. The method take 2-3hr per screen condition


%% some user definitions

% load the barcode to gene map
%load('/Users/mitchela2/Dropbox (UMass Medical School)/Mitchell Lab/Stocks-Bacteria/ASKA_Barcoded_Library/ASKA_lookup_map.mat');
load('./ASKA_lookup_map_LOCAL_COPY.mat'); % loacal copy, might not be updated

map = mapOdd; % the barcode lookup table. Use mapEven or mapOdd or map (for both)
sqCutoff = 10; % score quality cutoff
pattern_PS1 = 'AGCTGCTTCG'; % last 10 bp of PS1 (used as anchor to identify the barcode)
pattern_PS16 = 'GAATCTTCG'; % first 10 bp of PS16
%% iterate over files and extract barcode data on the fly
fastqList =dir('*.fastq');

allBarcodes = map.keys;

dataset = [];
parfor iFile=1:length(fastqList)
    tic;
    % Import data from the fastq file
    [Header, Sequence, Qual] = fastqread(fastqList(iFile).name);
    Header = []; % this variable is junk
    
    % initialize empty variables to hold the data
    n = length(Sequence);
    hits = {}; nPos = 0; nNeg = 0;
    hits.names = unique(ASKA.geneName);
    hits.counts = zeros(size(hits.names));
    hits.barcodeCounts = zeros(size(ASKA.barcode));
    for i=1:n
        curBarcode = '';
        % get the sequence and mask low quality region
        curSeq = Sequence{1,i};
        %kNUPC = strfind(curSeq,'CGCATGGATCGTACAAGTCG');
        %if(kNUPC)
        %    disp 'found nupC barcode';
        %    fprintf("%d %s\n", kNUPC, curSeq);
        %end
        if(length(curSeq)>=50)
            curQS_str = Qual{1,i};
            curQS = double(curQS_str) - 33; % covert string to interger quality score
            tf = (curQS<sqCutoff);
            curSeqMasked = curSeq;
            curSeqMasked(tf) = 'N';
            
            % METHOD (1): identify the barcode by PS1 anchor and extract 20bp
%             k = strfind(curSeq(1:30),pattern_PS1);  
%             curBarcode = seqrcomplement(curSeqMasked((k+10):(k+29)));
            
            % METHOD (2): identify barcode by looking for a match (all barcodes vs. region in curSeq)
            curBarcodeRegion = seqrcomplement(curSeqMasked(20:50)); % based on primer design the barcode should be at 25-45bp 
            tfValidBarcode = (cellfun(@length,allBarcodes)>15 & cellfun(@length,allBarcodes)<25); % barcode length must be between 15-25bp
            for iBarcodeQuery=1:length(map)
                if(tfValidBarcode(iBarcodeQuery))
                    curBarcodeQuery = allBarcodes{iBarcodeQuery};
                    k = strfind(curBarcodeRegion,curBarcodeQuery);
                    if(~isempty(k)) % the curBarcodeRegion matches the curBarcodeQuery
                        queryLength = length(curBarcodeQuery);
                        curBarcode = curBarcodeQuery;
                        break;
                    end
                end
            end
        else
            curBarcode = 'THE LENGTH OF THIS LINE IS TOO SHORT TO INCLUDE A BARCODE';
        end
        %disp(curBarcode)
        if(map.isKey(curBarcode))
            nPos = nPos+1;
            geneName = ASKA.geneName(map(curBarcode));
            hits.barcodeCounts(map(curBarcode)) = hits.barcodeCounts(map(curBarcode))+1;
            inx = find(strcmp(geneName,hits.names));
            hits.counts(inx) = hits.counts(inx)+1;
        else
            nNeg = nNeg+1;
        end
        if(~mod(i,10000))
            fprintf('File %d of %d / seq %d of %d\n',iFile,length(fastqList),i,n);
        end
    end
    % calculate the relative frequency of counts
    hits.freq = hits.counts./sum(hits.counts);
    hits.RPM = hits.freq*1000000;
    
    % save the hits variable into a array of structures called dataset
    dataset(iFile).fileName = fastqList(iFile).name;
    dataset(iFile).runtime = toc/60/60;
    dataset(iFile).hits = hits;
    dataset(iFile).nPositiveBarcodes = nPos;
    dataset(iFile).nNegativeBarcodes = nNeg;
    dataset(iFile).nReads = nPos+nNeg;
end

save dataset_QS10 dataset

