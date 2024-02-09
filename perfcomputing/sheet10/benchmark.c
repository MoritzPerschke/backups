#include "array.h"
#include "linked_list.h"
#include "unrolled_linked_list.h"
#include "tiered_array.h"
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// counters to count operations
int insertion_counter = 0;
int deletion_counter = 0;
int read_counter = 0;
int write_counter = 0;

// position in mix array
int pos_mix_array = 0;

// sum of the read operations in order to avoid compiler optimization
long sum;

void* thread_work(void* params) {
    pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
    struct_param* parameters = params;

    switch (parameters->data_structure) {
        case 0:
            while (1) {
                traverse_array((char*)parameters->ptr_datastructure,
                               parameters->numb_elements,
                               parameters->element_size, parameters->mix_array,
                               parameters->length_mix_array);
            }
        case 1:
            while (1) {
                traverse_linked_list(
                    (ll_node_t*)parameters->ptr_datastructure,
                    parameters->numb_elements, parameters->element_size,
                    parameters->mix_array, parameters->length_mix_array);
            }
        case 2:
            while (1) {
                traverse_unrolled_linked_list(
                    (ull_node_t*)parameters->ptr_datastructure,
                    parameters->numb_elements, parameters->element_size,
                    parameters->mix_array, parameters->length_mix_array);
            }
        case 3:
            while(1) {
                traverse_tiered_array(
                    (tiered_array_t*)parameters->ptr_datastructure,
                    parameters->numb_elements,
                    parameters->element_size,
                    parameters->mix_array,
                    parameters->length_mix_array
                );
            }
        default:
            return NULL;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: %s <data_structure> <instruction_mix> <element_size> "
               "<numb_elements>\n",
               argv[0]);
        exit(EXIT_FAILURE);
    }

    struct_param* parameters = malloc(sizeof(struct_param));

    // 0 for array, 1 for linked list, 2 for unrolled ll, 3 for tiered array
    parameters->data_structure = atoi(argv[1]);

    // values of ins/del
    parameters->instruction_mix = atoi(argv[2]);

    // element size in Bytes
    parameters->element_size = atoi(argv[3]);

    // number of elements
    parameters->numb_elements = atoi(argv[4]);

    // 0% ins/del 100% read/write
    Op mix_0[2] = { READ, WRITE };

    // 1% ins/del 99% read/write
    Op mix_1[200];

    int next_operation = 0;
    // initialize mix_1
    for (int i = 0; i < 200; i++) {
        if (i == 99) {
            mix_1[i] = 2;
        } else if (i == 199) {
            mix_1[i] = 3;
        } else {
            mix_1[i] = next_operation;
            next_operation = next_operation == 0 ? 1 : 0;
        }
    }

    // 10% ins/del 90% read/write
    Op mix_10[20];

    next_operation = 0;
    // initialize mix_2
    for (int i = 0; i < 20; i++) {
        if (i == 9) {
            mix_10[i] = 2;
        } else if (i == 19) {
            mix_10[i] = 3;
        } else {
            mix_10[i] = next_operation;
            next_operation = next_operation == 0 ? 1 : 0;
        }
    }

    // 50% ins/del 50% read/write
    Op mix_50[4] = { READ, INSERT, WRITE, DELETE };

    // to initialize pointer mix_array to choosen mix array and initialize
    // length_mix_array with correct length of the mix array
    switch (parameters->instruction_mix) {
        case 0:
            parameters->mix_array = mix_0;
            parameters->length_mix_array = 2;
            break;
        case 1:
            parameters->mix_array = mix_1;
            parameters->length_mix_array = 200;
            break;
        case 10:
            parameters->mix_array = mix_10;
            parameters->length_mix_array = 20;
            break;
        case 50:
            parameters->mix_array = mix_50;
            parameters->length_mix_array = 4;
            break;
        default:
            printf("unknown instruction mix %i\n", parameters->instruction_mix);
            exit(EXIT_FAILURE);
    }

    switch (parameters->data_structure) {
        case 0:
            parameters->ptr_datastructure =
                init_array(parameters->element_size, parameters->numb_elements);
            break;
        case 1:
            parameters->ptr_datastructure = init_linked_list(
                parameters->element_size, parameters->numb_elements);
            break;
        case 2:
            parameters->ptr_datastructure = init_unrolled_linked_list(
                parameters->element_size, parameters->numb_elements);
            break;
        case 3:
            parameters->ptr_datastructure = init_tiered_array(
                parameters->element_size,
                parameters->numb_elements
            );
        default:
            printf("unknown data_structure %i\n", parameters->data_structure);
            exit(EXIT_FAILURE);
    }

    // create worker thread and kill it after 10 seconds
    pthread_t tid_work;
    pthread_create(&tid_work, NULL, thread_work, parameters);
    sleep(10);
    pthread_cancel(tid_work);
    pthread_join(tid_work, NULL);

    switch (parameters->data_structure) {
        case 0:
            free_array(parameters->ptr_datastructure);
            break;
        case 1:
            free_linked_list(parameters->ptr_datastructure);
            break;
        case 2:
            free_unrolled_linked_list(parameters->ptr_datastructure);
            break;
        case 3:
            free_tiered_array(parameters->ptr_datastructure);
    }

    // sum is printed in order to prevent compiler optimation
    printf("sum: %ld\n", sum);

    printf("insertion_counter: %d\n", insertion_counter);
    printf("deletion_counter: %d\n", deletion_counter);
    printf("write_counter: %d\n", write_counter);
    printf("read_counter: %d\n", read_counter);

    printf("total operations:\n%d\n",
           read_counter + write_counter + insertion_counter + deletion_counter);
}