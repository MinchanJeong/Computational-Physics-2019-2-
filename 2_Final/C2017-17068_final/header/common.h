#ifndef COMMON_H
#define COMMON_H

struct allocd_list_struct{
	void *ptr;
	struct allocd_list_struct *prev;
};
typedef allocd_list_struct allocd_list;

extern allocd_list *chain;

extern int init_allocdlist();
extern int write_allocdlist(void*);
extern int free_allocdlist();

#endif
