// Using a begin statement, a thread is created for each statement. This is asynchronous

begin writeln("This is one thread: ");
begin writeln("This is another thread");
begin writeln("And another thread");
