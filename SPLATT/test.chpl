
// Get pointers different parts of 1D array
// and store them as rows in a 2D matrix of ptrs
/*var arr : [0..9] int;
var arr2 : [0..1] c_ptr(int);
arr = 1;

writeln(arr);
writeln("");

arr2[0] = c_ptrTo(arr[0]);
arr2[1] = c_ptrTo(arr[5]);

for i in 0..1 {
    for j in 0..4 {
        arr2[i][j] = -12;
    }
}

for i in 0..1 {
    for j in 0..4 {
        write(arr2[i][j], " ");
    }
    writeln("");
}

writeln("");
writeln(arr);

// Get ptr to 2D matrix and access it
*/

/*var mat : [0..2, 0..2] int;
mat = -1;

var ptrMat = c_ptrTo(mat);
for i in 0..2 {
    for j in 0..2 {
        ptrMat[(i*3)+j] = 1;
    }
}

writeln(mat);*/

// Given a ptr to 2D matrix, get a pointer to somewhere inside it
/*var mat : [0..2, 0..2] int;
mat = -1;

var ptrMat = c_ptrTo(mat);
for i in 0..2 {
    for j in 0..2 {
        ptrMat[(i*3)+j] = 1;
    }
}

var t = ptrMat+(3);
for i in 0..2 {
    t[i] = 2;
}

t = ptrMat+(3*2);
for i in 0..2 {
    t[i] = 3;
}

writeln(mat);*/

/*class foo {
    var buf : [0..3] int;
}

class bar {
    var fooArr : [0..1] foo;
}

var t = new bar();
for i in 0..1 {
    t.fooArr[i] = new foo();
    t.fooArr[i].buf = -1;
}

var x = t.fooArr[0];
x.buf = 100;

for i in 0..1 {
    writeln(t.fooArr[i].buf);
}*/


/*var buf : [0..2] c_ptr(real);

var t : [0..3] real;
t = 1.0;
var s : [0..3] real;
s = 2.0;
var u : [0..3] real;
u = 3.0;

buf[0] = c_ptrTo(t);
buf[1] = c_ptrTo(s);
buf[2] = c_ptrTo(u);

for i in 0..2 {
    for j in 0..3 {
        writeln(buf[i][j]);
    }
}*/

use LinearAlgebra;

var a = Matrix([1, 2, 3],
               [4, 5, 6],
               [7, 8, 9]);

writeln(a);
writeln("");
a /= 2;

writeln(a);
