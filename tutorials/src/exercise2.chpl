// Use reduce to compute the minimum and second minimum of the following array:

var D: domain(1) = (1..10);
var A: [D] int = (-5,6,-2012,-75,2012,48,-700,65,100,0);

var minVal: int = min reduce A;
var secondMinVal: int = min reduce ([i in D] if A[i] > minVal then A[i]);
writeln("The minimum of A is: ", minVal);
writeln("The second minimum of A is: ", secondMinVal);
