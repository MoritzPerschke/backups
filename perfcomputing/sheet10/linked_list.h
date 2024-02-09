#include "benchmark.h"

typedef struct ll_node {
    struct ll_node* prev;
    struct ll_node* next;
    // flexible array member (avoid additional allocation)
    char data[];
} ll_node_t;

ll_node_t* init_linked_list(int element_size, int numb_elements);

void free_linked_list(ll_node_t* head);

void traverse_linked_list(ll_node_t* head, int numb_elements, int element_size,
                          Op mix_array[], int length_mix_array);