    clear all, clc
    A = 1:10;
    B = 11:20;
    C = 21:30;
    X = [A; B; C];

    
    fileName = fopen('data/data.txt','w');
    %for j = 1:size(data,1)
        fprintf(fileName, 'File starts Here');
        fprintf(fileName, '%d,', A);
        %fprintf(fileName, '%d,', data(j,:)); 
        fprintf(fileName, '\n');
    %end
    fclose('all');