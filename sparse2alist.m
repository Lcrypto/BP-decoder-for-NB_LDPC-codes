function sparse2alist(Q, H, file_name)
    [M, N] = size(H);
    file = fopen(file_name, 'w');
    
    if Q==2
        fprintf(file, '%d %d\n', N, M);
    else
        fprintf(file, '%d %d %d\n', N, M, Q);
    end
    
    column_weights = sum(full((H>0)));
    row_weights = sum(full((H>0)'));
    c_max = max(column_weights);
    r_max = max(row_weights); 
    
    fprintf(file, '%d %d\n', c_max, r_max);
    
    fprintf(file, '%d ', column_weights(1:end-1));
    fprintf(file, '%d\n', column_weights(end));
    
    fprintf(file, '%d ', row_weights(1:end-1));
    fprintf(file, '%d\n', row_weights(end));
    
    for i = 1:N
        indexes = (find(H(:,i)))';
        values = (full(H(indexes, i)))';
        if Q == 2
            temp = indexes;
        else
            temp = [indexes; values];
        end
        fprintf(file, '%d ', temp);
        if length(indexes) < c_max
            if Q == 2
                fillers = zeros(1, c_max - length(indexes));
            else
                fillers = zeros(1, 2*(c_max - length(indexes)));
            end
            fprintf(file, '%d ', fillers);
        end
        fprintf(file, '\n');
    end
    
    for i = 1:M
        indexes = find(H(i,:));
        values = full(H(i, indexes));
        if Q == 2
            temp = indexes;
        else
            temp = [indexes; values];
        end
        fprintf(file, '%d ', temp);
        if length(indexes) < r_max
            if Q == 2
                fillers = zeros(1, r_max - length(indexes));
            else
                fillers = zeros(1, 2*(r_max - length(indexes)));
            end
            fprintf(file, '%d ', fillers);
        end
        fprintf(file, '\n');
    end
    
    fclose(file);
end
