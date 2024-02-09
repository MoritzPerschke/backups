#ifndef BENCHMARK_H
#define BENCHMARK_H

typedef enum { READ, WRITE, INSERT, DELETE } Op;

// variables to hold choosen parameters
typedef struct parameters {
    // choosen data structure
    int data_structure;

    // ratio between insertion/deletion/read/write
    int instruction_mix;

    // element size in Bytes
    int element_size;

    // number of elements stored in the array/linked list
    int numb_elements;

    // choosen mix array
    Op* mix_array;

    // length of mix_array
    int length_mix_array;

    // pointer to datastructure
    void* ptr_datastructure;
} struct_param;

// counters to count operations
extern int insertion_counter;
extern int deletion_counter;
extern int read_counter;
extern int write_counter;

// position in mix array
extern int pos_mix_array;

// sum of the read operations in order to avoid compiler optimization
extern long sum;

#endif