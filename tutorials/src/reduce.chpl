use Time;

var t: Timer;

const numRect = 10000000;
const D : domain(1) = 1..numRect;
const width = 2.0 / numRect; //rectangle width
const baseX = -1 - width/2; //baseX+width is midpoint of 1st rectangle

proc rectangleArea(i : int) { //computes area of rectangle i
    const x = baseX + i*width;
    return width * sqrt(1.0 - x*x);
}

var halfPI : real;

t.start();
halfPI = + reduce rectangleArea(D);
t.stop();

writeln("Result: ",2*halfPI);
writeln("Time: ", t.elapsed(), " seconds");
