%%% Compute R
%[distanceVector,R] = sampling(Exp_New_Ds_5k,6800);
%[falsePositiveVector,falseNegativeVector] = computeBestR(lsh,AQS_orig_plus_minus{1},AQS_dim_plus_minus{1});


%%% Initialize LSH
%lsh = initializeLSH(Exp_New_Ds_5k,size(Exp_New_Ds_5k,2),size(Exp_New_Ds_5k,1),R,10,55);

R=0.08;
lsh.R = R;

disp("Query Set");
usedQerySet = AQS_orig;

fp = zeros(1,size(usedQerySet,2));
neighboursInRDistance = {};

for i=1:17
    %%% params: Query Set, LSH obj, Number of querys, Break on first neighbor in R distance
    [neighboursInRDistance,neighbours,fp(i),result] = query(usedQerySet{i},lsh,size(usedQerySet{i},2),1);
end

%%% Compute Parameters L and k
%[k,L] = computeParameters(Exp_New_Ds,Exp_New_Qs,0.3,0.1,4);

function lsh = initializeLSH(dataSet,tableSize,dimensions,R,k,L)

    %L = 20;                     % Number of G hash functions and buckets
    M = dimensions;             % D-dimensions of the a element in dataset
    W = 4;                      % Bucket width
    %k = 20;                      % number of hash functions
    r1 = randi([1,2^5],1,8);    % random integer used in t1 
    r2 = randi([1,2^5],1,8);    % random integer used in t2 

    lsh = Lsh(dataSet, L, M, W, k, r1, r2,R);

    As = readmatrix("./txt/As.txt");
    Bs = readmatrix("./txt/Bs.txt");
    R1s = readmatrix("./txt/R1s.txt");
    R2s = readmatrix("./txt/R2s.txt");
    as = As(1:k*L,:);
    bs = Bs(1:k*L,:);
    r1s = R1s(k,1:k);
    r2s = R2s(k,1:k);

    m=1;
    for i=1:L
        h={};
        index= 1;
        for j=m:k + m - 1
            h{index} = {as(j,:),bs(j,:)};
            index = index + 1;
        end
        m=m+1;
        g{i} = h;
    end

    lsh.g = g;
    lsh.r1 = r1s;
    lsh.r2 = r2s;
    lsh.tableSize = tableSize;
    lsh.dataset = dataSet;
    lsh = lsh.initialize(dataSet);
end
function [neighboursInRDistance,neighbours,fp,result] = query(queryDataSet,lsh,numberOfQuerys,breakVar)
    result = {}; 
    fp = 0;
    neighboursInRDistance = {};
    %t1s = [];
    timeElapsed = [];
    neighbours = {};
    
    for q=1:numberOfQuerys%size(queryDataSet,2)
        row = queryDataSet(:,q);
        
        %%% Measure query time, used in k parameter computation
        tic
        
        [neighboursInRDistance{q},neighbours{q},t1,ok]  = lsh.queryAnom(row,breakVar); 
        timeElapsed(end+1) = toc;
        
        %disp(["Query point", q , ", Query Time: ",timeElapsed(end)]);
        %t1s{1,q} = t1;
        
        if (ok == 1)
           result{end+1,2} = 'true';
           result{end,3} = q;
           fp = fp + 1;
        else
           result{end+1,2} = 'false';
           result{end,3} = q;
        end
     end
end

function [col,dataSet] = createDataSet(csvFile)
        table = readtable(csvFile);
        dataSet = table2array(table)';
        %dataset = unique(dataset', 'rows')';
        [~,col] = size(dataSet);
        dataSet = dataSet * 10;
        minMatrix = min(dataSet(:));
        dataSet = dataSet + abs(minMatrix) + 10;
end
function result = bfsearch(querySet, dataSet,nrOfItems)
    result = [];
    closestPoint = 0;
    for i=1:nrOfItems
        min = euclideanDistance(querySet(:,i),dataSet(:,1));
        for j=1:size(dataSet,2)
            a = querySet(:,i);
            b = dataSet(:,j);        
            ed = euclideanDistance(a,b);
            if (ed < min)
                min = ed;
                closestPoint = j;
            end
        end
        result(i,1) = min;
        result(i,2) = closestPoint;
    end
end
function result = euclideanDistance(pointA,pointB)
            result = 0;
            for i=1:size(pointA,2)
                result = result + (pointA(i) - pointB(i))^2;
            end
            result = sqrt(result);          
end  
function [distanceVector,output] = sampling(dataSet,nrOfSamples)
    distVector = [];
    finalDistVector = [];
    r = randi([1 size(dataSet,2)],1,nrOfSamples);
    
    for i=1:nrOfSamples-1
        for j=i+1:nrOfSamples
            distVector(end+1) = euclideanDistance(dataSet(:,r(i)),dataSet(:,r(j)));
        end
        
        distVector = sort(distVector);
        distVector = distVector(distVector ~= 0);
        distVector = distVector(1);
        
        finalDistVector = [finalDistVector distVector];
        distVector = [];
   end
    
    output = mean(finalDistVector);
    distanceVector = finalDistVector;
end
function output = normalizeInput(dataSet)
    %minX = [0,0,0,-100,0,-40,0,0,0,0,0,-40,0,0,0];
    minX = [0,0,0,-100,0,-40,0,0,0,0,0,-40,0,0,0];
    
    maxX = [10000,100,255,100,100,215,100,100,100,1,16383.75,215,511.75,511.75,511.75];
    
%     minX = [0    ,0  ,0  ,0  , 0 ,0   ,0  ,-40 ,0   ,0  ,0  ,0     ,0     ,0     ,0     ,0        ,0];
%     maxX = [10000,100,100,255,100,6000,100,215 ,100 ,100,100,511.75,511.75,511.75,511.75,16383.75, 200];
    for i=1:size(dataSet,1)
        row = dataSet(i,:);
        dataSet(i,:) = norm(row,minX(i),maxX(i));
    end
    
    output = dataSet;
end
function output = norm(row,minX,maxX)
    output = [];
    for i=1:size(row,2)
       output(i) =  (row(i) - minX) / (maxX - minX);
    end
end
function [k,L] = computeParameters(dataSet,querySet,R,err,w)
    %L = computeL(w,pi,x); 
   
    %k = 21;
    k = computeK(dataSet,querySet);
    L = computeL(R,w,k,err);
end
function L = computeL(R,w,k,err)
  
    %PDF = @(x) exp(-x.^2/2) / sqrt(2*pi);
    %P1 =@(t) integral(PDF(t/R) / R * (1-t/w), t, 0, w);
    
%     PDF = @(x) exp(-x.^2/2) ./ sqrt(2*pi);
%     fun = @(t) (PDF(t/R)./R).*(1-t./w);


    %CDF = @(x) exp(-x.^2/2) ./ sqrt(2*pi);
    %normp = integral(CDF, -inf, w/R);
    normp1 = normcdf(-w);
    normp2 = normcdf(-w/c);
    P1 = 1 - 2*normp1 - (2/(sqrt(2*pi) * w))* (1 - exp(-(w^2/2)));
    P2 = 1 - 2*normp2 - (2/(sqrt(2*pi) * w/c))* (1 - exp(-(w^2/(2*c^2))));
    L=ceil(log(1/err)/(-log(1-P1^k)));
    
end
function p = computeP(w,c)
    P1 = 1 - 2*normp1 - (2/(sqrt(2*pi) * w))* (1 - exp(-(w^2/2)));
    P2 = 1 - 2*normp2 - (2/(sqrt(2*pi) * w/c))* (1 - exp(-(w^2/(2*c^2))));
    p = log(P1) / log(P2);
end
function k = computeK(baseDataSet,samSet)
    %Tg = compute L functions And Retrieve Buckets
    %Tc = time for computing distances to all points in the bucket
    %Return the k for which Tg + Mean(Tc) is minimal.
    Tg = [];
    Tc = [];
    TcMeanVector = [];
    SumVector = [];
    MinimumVector = [];
    L = 5;
    k = 0;
    while (k < 50)
        k=k+1;
        lsh = initializeLSH(baseDataSet,size(baseDataSet,2),size(baseDataSet,1),0,k,L);
        
        for i=1:size(samSet,2)
            row = samSet(:,i);
            %get Tc and Tg vectors for Q point
            [Tg(i,k),Tc(i,k)]  = lsh.contructLFunctionsAndComputeDistances(row);
            
            %Compute mean of Tc
            TcMeanVector(i) = mean(Tc(i,k));  
            SumVector(i,k) = Tg(i,k) + TcMeanVector(i);

            [~,I] = min(SumVector(i,:));
            MinimumVector(i) = I;
        end
    end
    k = floor(mean(MinimumVector));
end
function w = computeOptimalW()
     c=10;
     pVector = [];
     for w=1:10
         normP1 = normcdf(-w);
         normP2 = normcdf(-w/c);
         P1 = 1 - 2*normP1 - (2/(sqrt(2*pi) * w))* (1 - exp(-(w^2/2)));
         P2 = 1 - 2*normP2 - (2/(sqrt(2*pi) * w/c))* (1 - exp(-(w^2/2*c^2)));
         p = log(1/P1) / log(1/P2);
         pVector(end+1,1) = w;
         pVector(end,2) = p;
         pVector(end,3) = P1;
         pVector(end,4) = P2;
     end
     
end
function subset = getSubSet(dataSet,nrOfSamples)
    subset = [];
     r = randi([1 size(dataSet,2)],1,nrOfSamples);
     for i=1:nrOfSamples
         subset(:,end+1) = dataSet(:,r(i));
     end
end

function [falsePositiveVector,falseNegativeVector] = computeBestR(lsh,dataSet,anomalySet)
%%% Compute fp and fn rates for R ranging from 0.05 to 0.5 in steps of 0.01
    %%% vector form (R, value)
    falsePositiveVector = [];
    falseNegativeVector = [];
    
    for r=0.03:0.01:0.5
        disp(["Query: R=",r]); 
        lsh.R = r;
        falsePositiveVector(end+1,1) = r;
        falseNegativeVector(end+1,1) = r;
        [neighboursInRDistance,neighbours,falsePositiveVector(end,2),result] = query(dataSet,lsh,size(dataSet,2),1);
        [neighboursInRDistance,neighbours,falseNegativeVector(end,2),result] = query(anomalySet,lsh,size(anomalySet,2),1);
    end
    falsePositiveVector(:,2) = 500 - falsePositiveVector(:,2);
    falsePositiveVector(:,2) = falsePositiveVector(:,2) * 100 / 500;
    
    falseNegativeVector(:,2) = falseNegativeVector(:,2) * 100 / 500;
    
end