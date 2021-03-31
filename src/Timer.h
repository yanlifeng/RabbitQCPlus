#ifndef TIMER
#define TIMER
#include<time.h>
#include<sys/timeb.h>
#include<stdio.h>
#define TPRINT(start_t,end_t,str)\
    printf ("%s time: %6.3fs \n", \
	        str,\
            (end_t.time + end_t.millitm / 1000.0) - (start_t.time + start_t.millitm / 1000.0));
#endif