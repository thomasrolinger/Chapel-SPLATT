// No code runs until the cobegin block is finished

cobegin {
    writeln("This is one thread: ");
    writeln("This is another thread");
    writeln("And another thread");
}
begin writeln("All other threads are done");
