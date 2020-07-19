#include <stdlib.h>
#include <stdio.h>

#include "common.h"

allocd_list *chain;

int init_allocdlist(){
	
	chain = (allocd_list*)malloc(sizeof(allocd_list));

	chain->ptr =NULL;
	chain->prev = NULL;
	return 0;//
}

int write_allocdlist(void* p){

	allocd_list *_c = (allocd_list*)malloc(sizeof(allocd_list));
	_c->ptr = p;
	_c->prev = chain;
	chain = _c;
	return 0;
}

int free_allocdlist(){

	fprintf(stdout,"\nfreeing unfreed allocated memory");

	int count = 0;

	allocd_list *pivot;
	while(chain->prev != NULL){
		free(chain->ptr);
		pivot = chain;
		chain = chain->prev;
		free(pivot);
		count++;
	}

	fprintf(stdout,"\n");
	free(chain);
	fprintf(stdout,"%d Allocated pointers freed properly\n",count);
	return 0;
}
