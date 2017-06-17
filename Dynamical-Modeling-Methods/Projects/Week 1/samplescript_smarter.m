data = dlmread('sampledata1.DAT') ; 
time = data(:,1) ;
[minimum,index] = min(abs(time-50)) ;
f1 = data(:,2) ;
f1_norm = normalize(f1,index) ;
    
f2 = data(:,3) ;
f2_norm = normalize(f2,index) ;

f3 = data(:,4) ;
f3_norm = normalize(f3,index) ;

figure
hold on 
plot(time,f1_norm,'b')
plot(time,f2_norm,'g')
plot(time,f3_norm,'r')
