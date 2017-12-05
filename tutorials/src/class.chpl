// Classes are like pointers to structs. Assignment to a class copies the
// reference to the class object

class circle {
    var r: real;
    proc area() {
        return 3.14 * (r**2);
    }
}

var c1, c2: circle;
c1 = new circle(12.0);
c2 = c1; // makes another references to c1
c2.r = 10; // also changes c1.r

writeln("The area of c1 is ", c1.area());
writeln("The area of c2 is ", c2.area());

delete c1; // Need to free memory
