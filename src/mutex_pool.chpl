/*
    Author: Thomas Rolinger (tbrolin@cs.umd.edu)
    Date:   12/20/2017
    File:   mutex_pool.chpl

    Description:    This is a module file for things related
                    to mutexes
*/

module MutexPool {
    use Base;

    /*****************************
    *
    *   mutex_pool class
    *
    ******************************/
    /*
        Chapel doesn't have a literal mutex/lock
        mechanism, as far as I know. However, it
        has sync variables, which we can use as
        locks. When you attempt to write to a
        sync variable, you are getting the lock.
        And then when you read from the sync variable,
        you are releasing the lock.

        So to implement a mutex pool, we will have
        an array of sync variables. Also, like SPLATT,
        we will have padding between each sync variable
        to avoid false sharing. False sharing is when
        threads modify variables that reside on the same
        cache line. SPLATT stores locks as 32-bit ints but
        we'll use 64-bit ints. Therefore, to allow for 64-bytes
        of padding between each lock, our pad size is 8 bytes.
        So lock ID 0 is at locks[0] and lock ID 1 is at locks[8].
        The number of bytes between locks[0] and locks[8] is
        sizeof(int64)* 8 = 64 bytes.
    */

    class mutex_pool {
        var num_locks : int;
        var pad_size : int;
        var locks_d : domain(1) = 0..(DEFAULT_NLOCKS * DEFAULT_LOCK_PAD)-1;
        //var locks : [locks_d] sync bool;
        var locks : [locks_d] atomic bool;
    }

    /*****************************
    *
    *   Public Functions
    *
    ******************************/
    
    proc mutex_alloc(const num_locks : int, const pad_size : int)
    {
        var pool : mutex_pool = new mutex_pool();
        pool.num_locks = num_locks;
        pool.pad_size = pad_size;
        for l in 0..num_locks-1 {    
            const lock = mutex_translate_id(l, num_locks, pad_size);
            //pool.locks[lock].writeXF(true);
            pool.locks[lock].clear();
        }
        return pool;
    }

    /*########################################################################
    #   Description:    Convert an arbitrary integer ID to a lock ID in a
    #                   mutex pool.
    #
    #   Parameters:     id (int):           An abitrary integer ID (i.e. matrix row).
    #                   num_locks (int):    Size of mutex pool
    #                   pad_size (int):     Padding between each lock
    #
    #   Return:         int : lock ID
    ########################################################################*/
    proc mutex_translate_id(const id : int, const num_locks : int, const pad_size : int) : int
    {
        return (id % num_locks) * pad_size;
    }

    /*########################################################################
    #   Description:    Claim a lock of a mutex pool. The lock is identified
    #                   with an ID, which is an arbitrary integer which 
    #                   uniquely identifies the memory to protect. The ID
    #                   is then translated into an actual lock based on
    #                   the pool size and padding.
    #
    #   Parameters:     pool (mutex_pool):  The pool to use
    #                   id (int):           Lock ID
    #
    #   Return:         None
    ########################################################################*/
    proc mutex_set_lock(const pool : mutex_pool, const id : int, tid)
    {
        // Setting a lock in our case means reading the sync variable. This
        // sets the sync var to empty and means no other task will be able
        // to "set the lock" until I unset it (i.e. write to it).
        const lock_id = mutex_translate_id(id, pool.num_locks, pool.pad_size);
        //pool.locks[lock_id];
        while pool.locks[lock_id].testAndSet() {
            chpl_task_yield();
        }
    }

    /*########################################################################
    #   Description:    Release a lock set with mutex_set_lock()
    #
    #   Parameters:     pool (mutex_pool):  The pool containing the lock
    #                   id (int):           Lock ID
    #
    #   Return:         None
    ########################################################################*/
    proc mutex_unset_lock(const pool : mutex_pool, const id : int, tid)
    {
        // To unset the lock, we write to it. This sets it from empty to full.
        const lock_id = mutex_translate_id(id, pool.num_locks, pool.pad_size);
        //pool.locks[lock_id].writeXF(true);
        pool.locks[lock_id].clear();
    }
}
