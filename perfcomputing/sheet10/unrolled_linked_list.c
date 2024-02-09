#include "unrolled_linked_list.h"
#include "benchmark.h"
#include <stdlib.h>

ull_node_t* ull_create_node(int element_size) {
    ull_node_t* new_node = malloc(sizeof(ull_node_t) +
                                  element_size * sizeof(char) * ULL_ARRAY_LEN);

    // insert random values
    for (int i = 0; i < element_size * ULL_ARRAY_LEN; i++) {
        new_node->data[i] = i % 255;
    }
    new_node->elements = ULL_ARRAY_LEN;
    new_node->next = NULL;
    new_node->prev = NULL;
    return new_node;
}

void ull_free_node(ull_node_t* node) {
    free(node);
}

ull_node_t* init_unrolled_linked_list(int element_size, int numb_elements) {
    // create head of the linked_list
    ull_node_t* new_node = ull_create_node(element_size);
    ull_node_t* head = new_node;
    ull_node_t* tail = new_node;

    // create rest of the nodes
    for (int i = 1; i < numb_elements / ULL_ARRAY_LEN; i++) {
        ull_node_t* new_node = ull_create_node(element_size);
        new_node->prev = tail;
        tail->next = new_node;
        tail = new_node;
    }

    int elements = numb_elements % ULL_ARRAY_LEN;
    if (elements != 0) {
        ull_node_t* new_node = ull_create_node(element_size);
        new_node->prev = tail;
        tail->next = new_node;
        tail = new_node;
        // set element count to number of needed elements
        new_node->elements = elements;
    }

    return head;
}

void free_unrolled_linked_list(ull_node_t* head) {
    ull_node_t* current = head;
    ull_node_t* temp;

    while (current != NULL) {
        temp = current;
        current = current->next;
        ull_free_node(temp);
    }
}

void ull_read(ull_node_t* node, int index, int element_size) {
    for (int i = index * element_size; i < (index + 1) * element_size; i++) {
        sum += node->data[i];
    }
    read_counter++;
}

void ull_write(ull_node_t* node, int index, int element_size) {
    for (int i = index; i < (index + 1) * element_size; i++) {
        node->data[i] = i % 255;
    }
    write_counter++;
}

void ull_insert(ull_node_t* current, int index, int element_size) {
    if (current->elements < ULL_ARRAY_LEN) {
        // insert into array
        for (int i = element_size * (current->elements + 1) - 1;
             i >= element_size * index; i--) {
            current->data[i] = current->data[i - element_size];
        }

        for (int i = index * element_size; i < (index + 1) + element_size;
             i++) {
            current->data[i] = i % 255;
        }

        current->elements += 1;

        insertion_counter++;
    } else {
        // insert new node
        ull_node_t* new_node = ull_create_node(element_size);
        new_node->next = current->next;
        new_node->prev = current;

        ull_node_t* temp = current->next;
        if (temp != NULL) {
            temp->prev = new_node;
        }
        current->next = new_node;

        new_node->elements = 1;

        insertion_counter++;
    }
}

void ull_delete(ull_node_t* current, int index, int element_size) {
    if (current->elements > 1) {
        // remove element from array
        for (int i = index * element_size;
             i < (ULL_ARRAY_LEN - 1) * element_size; i++) {
            current->data[i] = current->data[i + element_size];
        }

        current->elements -= 1;

        deletion_counter++;
    } else {
        // remove node
        ull_node_t* previous = current->prev;
        ull_node_t* next = current->next;
        next->prev = previous;
        previous->next = next;

        free(current);
        deletion_counter++;
    }
}

void traverse_unrolled_linked_list(ull_node_t* head, int numb_elements,
                                   int element_size, Op mix_array[],
                                   int length_mix_array) {
    for (int i = 0; i < numb_elements; i++) {
        // to find current node
        ull_node_t* current = head;
        int prev_elements = 0;
        while (prev_elements + current->elements <= i) {
            prev_elements += current->elements;
            current = current->next;
        }

        int index = i - prev_elements;

        switch (mix_array[pos_mix_array]) {
            case READ:
                ull_read(current, index, element_size);
                break;
            case WRITE:
                ull_write(current, index, element_size);
                break;
            case INSERT:
                ull_insert(current, index, element_size);
                break;
            case DELETE:
                ull_delete(current, index, element_size);
                break;
        }
        pos_mix_array = (pos_mix_array + 1) % length_mix_array;
    }
}