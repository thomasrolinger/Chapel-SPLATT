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
        Mutex pool is implemented as an array of atomic boolean
        variables. We tried sync variables but with qthreads and our
        very short critical sections, they are too heavy-weight.
    */

    class mutex_pool {
        var num_locks : int;
        var pad_size : int;
        var locks_d : domain(1) = 0..(DEFAULT_NLOCKS * DEFAULT_LOCK_PAD)-1;
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
        const lock_id = mutex_translate_id(id, pool.num_locks, pool.pad_size);
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
        pool.locks[lock_id].clear();
    }
}
