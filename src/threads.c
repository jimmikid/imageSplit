//
//  threads.c
//  TesiProject
//
//  Created by Gabriele on 17/01/16.
//  Copyright © 2016 Gianmarco Stinchi. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <threads.h>


#if defined( __APPLE__ ) || defined( __linux )

        #define __DEF_THREAD_POSIX
        #include <memory.h>
        #include <pthread.h>

#ifndef EDEADLK
    #define EDEADLK         35      /* Resource deadlock would occur */
#endif

#ifndef EINVAL
	#define EINVAL          22      /* Invalid argument */
#endif

#ifndef ESRCH	
    #define ESRCH            3      /* No such process */
#endif


        struct thread
        {
            pthread_t m_thread;
        };

        thread* execute_task(task function,void* context)
        {
            //output
            thread* th =  (thread*)malloc(sizeof(struct thread));
            //...
            if(pthread_create(&th->m_thread,
                              NULL,
                              function,
                              context)  != 0)
            {
                free(th);
                return NULL;
            }
            
            return th;
        }

        void* joint(thread* th)
        {
            void* return_value = NULL;
            
            //get error
            int error = pthread_join(th->m_thread, &return_value);
            
            
            
            if (error != 0)
            {
                switch (error)
                {
                    case EDEADLK:  printf("EDEADLK\n"); break;
                    case EINVAL:   printf("EINVAL\n");  break;
                    case ESRCH:    printf("ESRCH\n");   break;
                    default: break;
                }
                return NULL;
            }
            
            return return_value;
        }
#else 

 	#define __DEF_THREAD_WIN32
    #include <windows.h>
 	#include <process.h>

	typedef struct args
	{
		thread* m_thread;
		void*   m_context;
		task    m_task;
	}
	args;

	unsigned WINAPI __wrap_function(void* ptr_args)
	{
		args* v_args = (args*)ptr_args;
		bool ret = v_args->m_task(v_args->m_context);
		_endthreadex(ret ? 1 : 0);
		return 0;
	}

	struct thread
	{
		HANDLE m_thread;
		args*  m_args;
	};

	thread* execute_task(task function, void* context)
	{
		//output
		thread* th = (thread*)malloc(sizeof(struct thread));
		//alloc args
		th->m_args = (args*)malloc(sizeof(struct args));
		//init args
		th->m_args->m_thread = th;
		th->m_args->m_context = context;
		th->m_args->m_task = function;
		//...
		if ( !(th->m_thread=(HANDLE)_beginthreadex(NULL, 0, __wrap_function,  th->m_args,  0, NULL)) )
		{
			free(th);
			return NULL;
		}

		return th;
	}

	void* get_task_context(thread* th)
	{
		return th->m_args->m_context;
	}

	bool joint(thread* th)
	{
		if (WaitForSingleObject(th->m_thread, INFINITE) == WAIT_FAILED) 
		{
			return false;
		}
		//get return value
		DWORD ret = 0;
		GetExitCodeThread(th->m_thread, &ret);
		//close thread
		CloseHandle(th->m_thread);
		//free structure
		free(th->m_args);
		free(th);

		return ret != 0;
	}
#endif


