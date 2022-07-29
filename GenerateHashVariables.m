maxK = 50;  % Maximum K H_i hash functions for each Gj
maxL = 100; % Maximum L G_j functions
M = 17;     % D-dimensions
W = 4;      % Bucket width

As = {};
Bs = {};
R1s = zeros(maxK);
R2s = zeros(maxK);

for index=1:maxK*maxL
    %As{index, 1} = 1 + 1000*randn(1, M);
    number = normrnd(0,4,1);
    x=[];
    for i=1:M
        number  = normrnd(0,4,1);
        while number < -4 || number > 4
           number  = normrnd(0,4,1);
        end
        x(end+1) = number;
    end
    As{index, 1} = x;
    Bs{index, 1} = unifrnd(0, W);
    
    %sqrr(-2^log(x1)*cos(2.0*pi*x2));
end

for i=1:maxK
    r1 =  randi([1,2^10],1,i);
    r2 = randi([1,2^10],1,i);
    [~,nrCol] = size(r1);
    
    for j=1:nrCol
        R1s(i,j) = r1(1,j);
        R2s(i,j) = r2(1,j);
    end
end

writematrix(cell2mat(As),"./txt/As.txt");
writematrix(cell2mat(Bs),"./txt/Bs.txt");
writematrix(R1s,"./txt/R1s.txt");
writematrix(R2s,"./txt/R2s.txt");