classdef Lsh
    %Implementation of Locality Sensitive Hashing Algorithm
  
    properties
        datasetFile  % Csv file containing the dataset
        dataset;     % Matrix of dataset. [D X N] N columns of D-dimensions
        L;           % Number of G hash functions and buckets
        M;           % D-dimensions of the a element in dataset
        D;           % Size of dataset
        tableSize;   % Data set Size
        W;           % Bucket width
        hashTables = {}; % Cell storing hash tables
        r1,r2;       % random numbers used in t1 and t2
        P;           % prime number used in t1 and t2
        k;           % number of hash functions
        g;           % g functions
        unifRVector;
        hashF = {};
        t1s = [];    % the t1 values for the g functions. (the kays)
        R;           % Max distance from query point
    end
    
    methods
        function obj = Lsh(datasetFile, L, M, W, k,r1,r2,R)
            obj.datasetFile = datasetFile;
            obj.L = L;
            obj.M = M;
            obj.W = W;
            obj.P = 2^32 - 5;
            obj.k = k;
            obj.g = {};
            obj.r1 = r1;
            obj.r2 = r2;
            obj.R = R;
        end   
        function lsh = initialize(this,dataset)

            %this.dataset = Lsh.parseCsv(this.datasetFile);
            %this.dataset = unique(this.dataset', 'rows')';
            this.dataset = dataset;
            this.D = size(this.dataset, 2);
            %this.tableSize = this.D;
            [this.t1s,this.hashF,this.hashTables] = this.generateBuckets();
            lsh = this;
            
        end
        function lsh = initializeWithDimensionIncrease(this)
            this.D = size(this.dataset, 2);
            [this.hashF,this.hashTables] = this.generateBuckets();
            lsh = this;
        end
        function out = generateGFunctions(this)
            out = {};
            for i=1:this.L
                 out{i} = generateHashFam(this); 
            end
        end
        function [neighbors, t1Vector,ok] = queryAut(this, query)
           ok = 1;
            output = 1;
            t1Vector = [];
            %used to check if all hash tables return true
            for i=1:this.L
                
                hashValue = [];
                
                for h=1:this.k
                    hashValue(h) = this.hashVector(query',this.g{i}{h});
                end
                
                t1 = this.calculateT1(hashValue); 
                t1Vector(i,1) = t1;

                if ~isKey(this.hashTables{i}, t1)
                    ok = 0;
                    neighbors{i} = 0;
                else
                    % get all points in the bucket
                    neighbors{i} = cell2mat(values(this.hashTables{i},{t1}));

                    % go over all neighbours
                    neihgborSize = size(neighbors{i},2);
                    
                    for idx = 1:neihgborSize
                        % compute euclidian distance from query to each
                        % point in the bucket
                        eudist = this.euclideanDistance(query',this.dataset(:,neighbors{i}(idx)));
                       
                        %check if distance from query to point is 0!
                        if ( eudist ~= 0)
                            ok = 0;
                        else
                            output = 0;
                        end
                    end 
                end
            end
            
            if (output == 0)
                ok = 1;
            end
            
        end
        function [neighborsInRDistance,neighbors, t1Vector,ok] = queryAnom(this, query, breakVar)
            ok = 1;
            neighborInRDistanceFound = 1;
            t1Vector = [];
            totalNeighboursInRDistance = 0;
            neighborsInRDistance = {};
            for i=1:this.L
                hashValue = [];
                
                for h=1:this.k
                    hashValue(h) = this.hashVector(query',this.g{i}{h});
                end
                
                t1 = this.calculateT1(hashValue); 
                t1Vector(i,1) = t1;
                               
                if ~isKey(this.hashTables{i}, t1)
                    ok = 0;
                    neighbors{i} = 0;
                    neighborsInRDistance{i} = 0;
                else
                    % get all points in the bucket
                    neighbors{i} = cell2mat(values(this.hashTables{i},{t1}));
                    
                    % go over all neighbours
                    neihgborSize = size(neighbors{i},2);
                    
                    for idx = 1:neihgborSize
                        % compute euclidian distance from query to each
                        % point in the bucket
                        eudist = this.euclideanDistance(query',this.dataset(:,neighbors{i}(idx)));
                        
                        %compare distance to R
                        if ( eudist > this.R)
                            ok = 0;
                           % disp(["Eudist: ",eudist, "- R: ",this.R]);
                        else
                            %%% at least one point in R distance exists
                            neighborInRDistanceFound = 0;
                            totalNeighboursInRDistance = totalNeighboursInRDistance + 1;
                             
                            if breakVar == 1
                                break;
                            end
                           
                        end
                    end 
                end
                neighborsInRDistance{i} = totalNeighboursInRDistance;
                totalNeighboursInRDistance = 0;
                
                if neighborInRDistanceFound == 0 && breakVar == 1
                    break;
                end
             
            end
            
            if (neighborInRDistanceFound == 0)
                ok = 1;
            else
                ok = 0;
            end
        end
        function [Tg,Tc] = contructLFunctionsAndComputeDistances(this,query)
            
            %compute all G functions
            tic
            for i=1:this.L
                hashValue = [];
                
                for h=1:this.k
                    hashValue(h) = this.hashVector(query',this.g{i}{h});
                end
                
                t1 = this.calculateT1(hashValue); 
                t1Vector(i,1) = t1;
                
                if isKey(this.hashTables{i}, t1)
                   % get all points in the bucket
                    neighbors{i} = cell2mat(values(this.hashTables{i},{t1}));
                else
                    neighbors{i} = cell2mat({0});
                end
            end
            Tg = toc;
            
            tic
            %Retrieve all points from buckets
            for i=1:this.L
                % go over all neighbours
                %disp( size(neighbors{i},2));
                neihgborSize = size(neighbors{i},2);
                
                if (neihgborSize ~= 1)
                    for idx = 1:neihgborSize
                        % compute euclidian distance from query to each
                        % point in the bucket
                        %disp(["i=",i," idx=",idx]);
                        eudist = this.euclideanDistance(query',this.dataset(:,neighbors{i}(idx)));   
                    end 
                end
            end
            Tc = toc;
        end
    end
    methods(Hidden) 
        function [t1s,hashF, hashTables] = generateBuckets(this)
            index = 1;
            hashTables = {};
            hashF = {};
            t1s = [];
            
            for i=1:this.L
                hashTable = containers.Map('KeyType','int32','ValueType','any');
                
                for j=1:this.D
                    hashValue = [];
                    for h=1:this.k
                        row = this.dataset(:,j);
                        hashValue(h) = this.hashVector(row', this.g{i}{h}); % build each g function
                    end
                                       
                    hashF{index} = hashValue;
                    index = index + 1;
                    
                    t1 = this.calculateT1(hashValue); %index in hash table
                    t1s(i,j) = t1;
                    
                    if not(isKey(hashTable, t1))
                        hashTable(t1) = j;
                    else
                        hashTable(t1) = [cell2mat(values(hashTable,{t1})) j];
                    end
                    
                end
            hashTables{i} = containers.Map();
            hashTables{i} = hashTable;
            end  
        end
        function [out] = hashVector(this,inputVector,g)
            a = g{1};
            b = g{2};
            out = floor((dot(a, inputVector) + b)/this.W);
        end
        function out = calculateT1(this, g)
            sum_ = sum(g .* this.r1);
            modP = mod(sum_, this.P);
            out = mod(modP, this.tableSize);
        end
        function out = generateHashFam(this)
            out = {};
            for i=1:this.k
                a = 100 + 30*randn(1,this.M);
                b = this.uniformRandomNumber();
                out{i} = {a,b};
            end
        end
        function out = uniformRandomNumber(this)
            r = randi([0 10000000]);
            out = this.unifRVector(r);            
        end
    end
    methods(Static)
        function dataset = parseCsv(csvFile)
            table = readtable(csvFile);
            dataset = table2array(table)';
        end
        function result = euclideanDistance(pointA,pointB)
            result = 0;
            for i=1:size(pointA,2)
                result = result + (pointA(i) - pointB(i))^2;
            end
            result = sqrt(result);          
        end  
    end
end

