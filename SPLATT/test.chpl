
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
var mat : [0..2, 0..2] int;
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

writeln(mat);
