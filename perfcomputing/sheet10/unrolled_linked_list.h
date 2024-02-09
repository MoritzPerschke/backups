#include "benchmark.h"
#include <stdlib.h>

#define ULL_ARRAY_LEN 8

typedef struct ull_node {
    struct ull_node* prev;
    struct ull_node* next;
    int elements;
    // flexible array member (avoid additional allocation)
    char data[];
} ull_node_t;

ull_node_t* init_unrolled_linked_list(int element_size, int numb_elements);

void free_unrolled_linked_list(ull_node_t* head);

void traverse_unrolled_linked_list(ull_node_t* head, int numb_elements,
                                   int element_size, Op mix_array[],
                                   int length_mix_array);