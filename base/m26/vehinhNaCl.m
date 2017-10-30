clear all
fid = fopen('e:/testbil/bo64E.p1');
s = fgets(fid)
ncol = fscanf(fid,'%*c %d') 
for i = 1:2
   s = fgets(fid)
end
ncol = ncol +1
p1 = fscanf(fid,'%e %e',[ncol, inf]); % It has two rows now.
p1 =p1';
fclose(fid)

n = size(p1,1)

fid = fopen('e:/testbil/NaCl','w');
for i = 1:500: n
    fprintf(fid,'%e %e\n',p1(i,1), p1(i,9));
end 
fclose(fid)

