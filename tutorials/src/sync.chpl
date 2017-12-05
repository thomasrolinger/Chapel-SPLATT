proc statement1() {
    begin writeln("I'm one thread");
    begin writeln("I'm another thread");
}

proc statement2() {
    begin writeln("I'm yet another thread");
}

// Since we call these functions from within a sync block, the tasks spawned
// by the functions will finish before we continue
//sync {
//    begin statement1();
//    begin statement2();
//}
//begin writeln("I am less important and can wait");

// However, using a cobegin like this, the program will not be forced to
// wait for the tasks launched inside the functions to finish
cobegin {
    statement1();
    statement2();
}
begin writeln("Maybe I am just as important now.");
