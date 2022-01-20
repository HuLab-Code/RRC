function [halfTime, strtInt, endInt] = calcHalfDecay(risingInts,risingTimes)

strtInt = risingInts(1);
endInt = risingInts(end);
halfInt = strtInt + (endInt - strtInt)/2;
diffVals = (risingInts - halfInt);
negDiffVals = diffVals<0;
ind1 = find(negDiffVals,1,'first');
posDiffVals  = diffVals>0;
ind2 = find(posDiffVals,1,'last');
t1 = risingTimes(ind1);
t2 = risingTimes(ind2);
int1 = risingInts(ind1);
int2 = risingInts(ind2);
halfTime = t1 + (halfInt - int1)*((t2-t1)/(int2-int1));