// When a forall loop runs, it creates a thread for each core of the processor

var sum : int = 0;
forall i in 1..1000000 with (ref sum) {
    sum += i;
}
writeln(sum);
