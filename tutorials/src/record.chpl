// A record is like a struct in C, but not a pointer to a struct.
// In other words, they are value-based and an assignment of a record
// copies the record.

record circle {
    var r: real;
    proc area() {
        return 3.14 * (r**2);
    }
}

var c1, c2: circle; // declare two circles
c1 = new circle(12.0); // create/initilize c1
c2 = c1; // literally copies the data from c1 to c2
c2.r = 10; // access the record variable

writeln("The area of c1 is ", c1.area());
writeln("The area of c2 is ", c2.area());
