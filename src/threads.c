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
        #define EDEADLK         35      /* Resource deadlock would occur */
        #define EINVAL          22      /* Invalid argument */
        #define ESRCH            3      /* No such process */


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

#endif


