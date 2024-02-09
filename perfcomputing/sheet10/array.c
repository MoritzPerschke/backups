#include "benchmark.h"
#include <stdlib.h>

char* init_array(int element_size, int numb_elements) {
    char* array = calloc(element_size * (numb_elements + 1), sizeof(char));
    for (int i = 0; i < numb_elements * element_size; i++) {
        array[i] = i % 255;
    }
    return array;
}

void free_array(char* array) {
    free(array);
}

void read_array(char* array, int index, int element_size) {
    for (int i = index; i < index + element_size; i++) {
        sum += array[i];
    }
    read_counter++;
}

void write_array(char* array, int index, int element_size) {
    for (int i = index; i < index + element_size; i++) {
        array[i] = i % 255;
    }
    write_counter++;
}

void insertion_array(char* array, int index, int element_size,
                     int numb_elements) {
    // shift elements
    for (int i = (element_size * (numb_elements + 1)) - 1; i >= index; i--) {
        array[i] = array[i - element_size];
    }

    // insert random values
    for (int i = index; i < index + element_size; i++) {
        array[i] = i % 255;
    }
    insertion_counter++;
}

void deletion_array(char* array, int index, int element_size,
                    int numb_elements) {
    // shift elements
    for (int i = index + numb_elements; i < element_size * (numb_elements + 1);
         i++) {
        array[i - element_size] = array[i];
    }
    deletion_counter++;
}

void traverse_array(char* array, int numb_elements, int element_size,
                    Op mix_array[], int length_mix_array) {
    for (int i = 0; i < numb_elements * element_size; i += element_size) {
        switch (mix_array[pos_mix_array]) {
            case READ:
                read_array(array, i, element_size);
                break;
            case WRITE:
                write_array(array, i, element_size);
                break;
            case INSERT:
                insertion_array(array, i, element_size, numb_elements);
                break;
            case DELETE:
                deletion_array(array, i, element_size, numb_elements);
                break;
        }
        pos_mix_array = (pos_mix_array + 1) % length_mix_array;
    }
}