#include "benchmark.h"

char* init_array(int element_size, int numb_elements);

void free_array(char* array);

void traverse_array(char* array, int numb_elements, int element_size,
                    Op mix_array[], int length_mix_array);