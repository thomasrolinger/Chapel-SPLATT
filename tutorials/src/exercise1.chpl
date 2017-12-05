// Calculates powers of 2 up to 2^n and stores them in an array.
// Then prints the array. 
config var n: int = 10;

var A: [1..n] int;

for i in {1..n} {
    A[i] = 2**i;
}

for a in A {
    writeln(a);
}
