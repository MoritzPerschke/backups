#include "benchmark.h"
#include <stdlib.h>
#include <string.h>

#define _DEPTH 1     // how many tiers
#define _CONTAINER_SIZE 8 // how many elements per tier

typedef struct tiered_array{
    int num_containers;
    char*** containers;
}tiered_array_t;

tiered_array_t* init_tiered_array(int element_size,
                                  int num_elements);

void free_tiered_array(char* array);

void traverse_tiered_array(tiered_array_t* array,
                           int num_elements,
                           int element_size,
                           Op pos_mix_array[],
                           int num_op);