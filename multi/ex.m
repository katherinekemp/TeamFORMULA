A = magic(4);

fileID = fopen('myfile.txt','w');
nbytes = fprintf(fileID,'%5d %5d %5d %5d\n',A)