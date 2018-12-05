#ifndef __PAGESET_H__
#define __PAGESET_H__

/*pageset stack, record page numbers in a stack*/
struct pageset {
    int *data;       /*page number stack*/
    int npages;      /*num of pages in the satck*/
    int datasize;    /*stack size*/
};



typedef struct pageset *pageset_t;

pageset_t ps_new();
pageset_t ps_dup(pageset_t ps);
void ps_delete(pageset_t ps);
void ps_move(pageset_t from, pageset_t to);

void ps_clear(pageset_t ps);
void ps_push(pageset_t ps, int page);
int ps_pop(pageset_t ps);
int ps_size(pageset_t ps);
int ps_get(pageset_t ps, int i);
void ps_set(pageset_t ps, int i, int page);
void ps_replace(pageset_t ps, int from, int to);
void ps_remove(pageset_t ps, int page);
void ps_removeset(pageset_t ps, pageset_t set);
void ps_randomise(pageset_t ps);



#endif //__PAGESET_H__
