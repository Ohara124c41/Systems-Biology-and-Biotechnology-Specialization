function norm = normalize( ...
    vector,numpoints)

sum = 0 ;
for i=1:numpoints
  sum = sum + vector(i) ;  
end
average = sum/numpoints ;
norm = vector/average ;

return
