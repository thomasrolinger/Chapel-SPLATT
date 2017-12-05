var A: [1..3] int = (5,3,9); // 1D array, 3 elements 
var B: [1..3, 1..5] real; // 2D array, 3x5

A = 1; // set all elements to 1
A = B; // copies elements of B to A
A = A + 1; // adds 1 to each element of the array
A = A * 2; // multiplies each element by 2
