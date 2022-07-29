L = 10; % Number of g functions
k = 10; % Number of h functions
M = 3; % Data dimensionality
w = 4;

%%% Run once: generateAandB(), generateR1();
% Generate a and b (used in h computation)
AB = generateAandB(k, L, M, w);

% Generate r1 vector (used in t1 computation)
r1 = generateR1(k);

% Each column contains the t1 value for each g function
t1 = [];

% Input Vector
v = [1, 2, 3; 
    1, 2, 3;
    1, 2, 3];

for i = 1:size(v,2)
    for j = 1:L
        % Input vector, r1 vector, k hash func, a, b, w
        t1(i,j) = computeT1(v(:,i), r1, k, AB{j,1}, AB{j,2}, w);
    end
end

% Data point v, vector r1, number of hash functions k, a vector, b scalar,w
function t1 = computeT1(v, r1, k, a, b, w)
            h = [];
            % For a given point compute all h functions
            for j = 1:k
                h(end+1) = computeH(v, a{j,1}, b{j,1}, w);
            end
            % compute t1 value
            t1 = sum(h .* r1);
end

function hash = computeH(inputVector, a, b, w)
            hash = floor((dot(a, inputVector) + b)/w);
end

function AB = generateAandB(maxK, maxL, M, w)
    As = {};
    Bs = {};
    AB = {};
    
    minimumA = -4;
    maximumA = 4;
    
    for i=1:maxL
        for j=1:maxK
            number = normrnd(0,4,1);
            x=[];
            for ii=1:M
                number  = normrnd(0,4,1);
                while number < minimumA || number > maximumA
                   number  = normrnd(0,4,1);
                end
                x(end+1) = number;
            end
            As{j, 1} = x;
            Bs{j, 1} = unifrnd(0, w);
        end
        AB{i,1} = As;
        AB{i,2} = Bs;
    end
end

function r1 = generateR1(k)
    for i=1:k
        r1 =  randi([1,2^10],1,i);
    end
end
