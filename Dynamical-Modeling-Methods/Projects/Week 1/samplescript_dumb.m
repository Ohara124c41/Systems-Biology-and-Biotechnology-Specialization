data = dlmread('sampledata1.DAT') ; 
time = data(:,1) ;
f1 = data(:,2) ; 
[minimum,index] = min(abs(time-50)) ;
sum = 0 ;
for i=1:index
  sum = sum + f1(i) ;  
end
f1avg = sum/index ;
f1_norm = f1/f1avg ;
    
f2 = data(:,3) ;
sum = 0 ;
for i=1:index
  sum = sum + f2(i) ;  
end
f2avg = sum/index ;
f2_norm = f2/f2avg ;

f3 = data(:,4) ;
sum = 0 ;
for i=1:index
  sum = sum + f3(i) ;  
end
f3avg = sum/index ;
f3_norm = f3/f3avg ;

figure
hold on 
plot(time,f1_norm,'b')
plot(time,f2_norm,'g')
plot(time,f3_norm,'r')

