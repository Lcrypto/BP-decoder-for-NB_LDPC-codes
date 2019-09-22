
function peg2lad( pegFile, ladFile )

tic

%loop until full rank H_square is generated
while 1
    fid = fopen( pegFile, 'r' );
    n = fscanf(fid, '%d', 1);
    m = fscanf(fid, '%d', 1);
    maxDegree = fscanf(fid, '%d', 1);

    expansion = n./m;

    H_square = zeros(n, n);

    for i=1:m
        oneLocations = fscanf(fid, '%d', maxDegree)';
        oneLocations = oneLocations(find(oneLocations));
        degree = length(oneLocations);
        oneLocations = oneLocations(randperm(degree));

        for j=1:expansion
            locationSubset = oneLocations( round((j-1)./expansion.*degree)+1:round(j./expansion.*degree) );
            H_square( (i-1).*expansion+j, locationSubset ) = 1;
        end
    end

    fclose(fid);

    H_square = H_square(:, randperm(n));

    % check whether H_square is full rank in GF(2)
    A = int8(H_square);
    rank = 0;
    for col=1:n
        oneRows = find(A(:,col)==1);
        if length(oneRows)>0
            rank = rank+1;
            for row=oneRows(length(oneRows):-1:1)'
                A(row,col:n) = xor(A(row,col:n), A(oneRows(1),col:n));
            end
        end
    end
    
    toc

    if rank==n
        break;
    end
end
clear A

numCodes = 65;
period = 66;
txSeq = [33 16 24 8 28 12 20 4 30 14 22 6 26 10 18 2 31 15 23 7 27 11 19 3 29 13 21 5 25 9 17 1 32];
txSeq = reshape([txSeq+33; txSeq], 1, 66);
% txSeq is the order in which samples of accumulated syndrome are transmitted

fid = fopen(ladFile, 'w');

fprintf(fid, '%d ', [numCodes n sum(sum(H_square)) period]);

jc = [0 full(cumsum(sum(H_square)))];
fprintf(fid, '\n');
fprintf(fid, '%d ', jc);

for code=period-numCodes+1:period
    H = zeros(n/period*code, n);

    row = 1;
    txSubseq = mod(txSeq(1:code), period);
    for row_square=1:n
        H(row, :) = H(row, :) + H_square(row_square, :);
        row = row + ismember( mod(row_square, period), txSubseq );
    end

    fprintf(fid, '\n');
    fprintf(fid, '%d ', code);

    fprintf(fid, '\n');
    fprintf(fid, '%d ', sort(txSeq(1:code))-1);

    ir = mod(find(H)-1, n/period*code);
    fprintf(fid, '\n');
    fprintf(fid, '%d ', ir);

    code
    toc
    
    clear H
end

fclose(fid);