clear all
fid = fopen('e:/testbil/bo64D.p1');
s = fgets(fid)
ncol = fscanf(fid,'%*c %d') 
for i = 1:2
   s = fgets(fid)
end
ncol = ncol +1
p1 = fscanf(fid,'%e %e',[ncol, inf]); % It has two rows now.
p1 =p1';
fclose(fid)

fid = fopen('e:/testbil/bo64D.p2');
s = fgets(fid)
ncol = fscanf(fid,'%*c %d') 
for i = 1:2
   s = fgets(fid)
end
ncol = ncol +1
p2 = fscanf(fid,'%e %e',[ncol, inf]); % It has two rows now.
p2 =p2';
fclose(fid)

fid = fopen('e:/testbil/bo64D.p3');
s = fgets(fid)
ncol = fscanf(fid,'%*c %d') 
for i = 1:2
   s = fgets(fid)
end
ncol = ncol +1
p3 = fscanf(fid,'%e %e',[ncol, inf]); % It has two rows now.
p3 =p3';
fclose(fid)

n = size(p1,1)
wB = zeros(n,2);
wB(1,2) =0.972; wB(1,1) =0.;
wP = zeros(n,2);
wP(1,2) =1.382; wP(1,1) =0.;

for i = 2: n
   wP(i,1) = p1(i,1);
   wP(i,2) = wP(i-1,2)+2.8353e-4*(p1(i,1)-p1(i-1,1))*0.5*(p1(i,6)+p1(i,7)+p1(i-1,6)+p1(i-1,7))*1e3;
   wP(i,2) = wP(i,2) - 2.8353e-4*(p1(i,1)-p1(i-1,1))*0.5*(p3(i,6)+p3(i,7)+p3(i-1,6)+p3(i-1,7))*1e3;
   wP(i,2) = wP(i,2) - 2.8353e-4*(p1(i,1)-p1(i-1,1))*0.5*(p2(i,6)+p2(i,7)+p2(i-1,6)+p2(i-1,7))*1e3;  
   wB(i,1) = p1(i,1);
   wB(i,2) = wB(i-1,2)+2.8353e-4*(p1(i,1)-p1(i-1,1))*0.5*(p3(i,6)+p3(i,7)+p3(i-1,6)+p3(i-1,7))*1e3;
   wB(i,2) = wB(i,2) + 2.8353e-4*(p1(i,1)-p1(i-1,1))*0.5*(p2(i,6)+p2(i,7)+p2(i-1,6)+p2(i-1,7))*1e3;     
end 

plot(wP(:,1),wP(:,2))
hold on
plot(wB(:,1),wB(:,2))

fid = fopen('e:/testbil/waterB','w');
for i = 1: n
    fprintf(fid,'%e %e\n',wB(i,1), wB(i,2));
end 
fclose(fid)

fid = fopen('e:/testbil/waterP','w');
for i = 1: n
    fprintf(fid,'%e %e\n',wP(i,1), wP(i,2));
end 
fclose(fid)