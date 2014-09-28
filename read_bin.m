clear  

nx = 61;
ny = 61;
nz = 61;

f1 = fopen('F1.dat','rb');
f2 = fopen('F2.dat','rb');

F1 = fread(f1,nz*nx*ny,'real*4');
F2 = fread(f2,nz*nx*ny,'real*4');

F1 = reshape(F1,[nx ny nz]);
F2 = reshape(F2,[nx ny nz]);

fclose('all')