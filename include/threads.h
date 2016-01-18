//
//  threads.h
//  TesiProject
//
//  Created by Gabriele on 17/01/16.
//  Copyright © 2016 Gianmarco Stinchi. All rights reserved.
//

#ifndef threads_h
#define threads_h
//task
typedef void* (*task)(void*) ;

//thread
struct thread;
typedef struct thread thread;

thread* execute_task(task function,void* context);
void* joint(thread*);

#endif /* threads_h */
