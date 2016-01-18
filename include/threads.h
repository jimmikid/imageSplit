//
//  threads.h
//  TesiProject
//
//  Created by Gabriele on 17/01/16.
//  Copyright © 2016 Gianmarco Stinchi. All rights reserved.
//

#ifndef threads_h
#define threads_h
//boolean lib
#include <stdlib.h>
#include <stdbool.h>
//task
typedef bool (*task)(void*) ;

//thread
struct thread;
typedef struct thread thread;

thread* execute_task(task function,void* context);
void* get_task_context(thread*);
bool joint(thread*);

#endif /* threads_h */
