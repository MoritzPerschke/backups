#include "linked_list.h"
#include "benchmark.h"
#include <stdlib.h>

ll_node_t* ll_create_node(int element_size) {
    ll_node_t* new_node =
        malloc(sizeof(ll_node_t) + element_size * sizeof(char));

    // insert random values
    for (int i = 0; i < element_size; i++) {
        new_node->data[i] = i % 255;
    }
    new_node->next = NULL;
    new_node->prev = NULL;
    return new_node;
}

void ll_free_node(ll_node_t* node) {
    free(node);
}

ll_node_t* init_linked_list(int element_size, int numb_elements) {
    // create head of the linked_list
    ll_node_t* new_node = ll_create_node(element_size);
    ll_node_t* head = new_node;
    ll_node_t* tail = new_node;

    // create rest of the nodes
    for (int i = 1; i < numb_elements; i++) {
        ll_node_t* new_node = ll_create_node(element_size);
        new_node->prev = tail;
        tail->next = new_node;
        tail = new_node;
    }

    return head;
}

void free_linked_list(ll_node_t* head) {
    ll_node_t* current = head;
    ll_node_t* temp;

    while (current != NULL) {
        temp = current;
        current = current->next;
        ll_free_node(temp);
    }
}

void ll_read(ll_node_t* node, int element_size) {
    for (int i = 0; i < element_size; i++) {
        sum += node->data[i];
    }
    read_counter++;
}

void ll_write(ll_node_t* node, int element_size) {
    for (int i = 0; i < element_size; i++) {
        node->data[i] = i % 255;
    }
    write_counter++;
}

void ll_insert(ll_node_t* current, int element_size) {
    // insert new node after current node
    ll_node_t* new_node = ll_create_node(element_size);
    new_node->next = current->next;
    new_node->prev = current;

    // update pointer of old nodes
    ll_node_t* temp = current->next;
    if (temp != NULL) {
        temp->prev = new_node;
    }
    current->next = new_node;

    // insert random values
    for (int i = 0; i < element_size; i++) {
        new_node->data[i] = i % 255;
    }
    insertion_counter++;
}

void ll_delete(ll_node_t* current) {
    // update pointer
    ll_node_t* previous = current->prev;
    ll_node_t* next = current->next;
    next->prev = previous;
    previous->next = next;

    free(current);
    deletion_counter++;
}

void traverse_linked_list(ll_node_t* head, int numb_elements, int element_size,
                          Op mix_array[], int length_mix_array) {
    for (int i = 0; i < numb_elements; i++) {
        // to find current node
        ll_node_t* current = head;
        for (int j = 0; j < i; j++) {
            current = current->next;
        }

        switch (mix_array[pos_mix_array]) {
            case READ:
                ll_read(current, element_size);
                break;
            case WRITE:
                ll_write(current, element_size);
                break;
            case INSERT:
                ll_insert(current, element_size);
                break;
            case DELETE:
                ll_delete(current);
                break;
        }
        pos_mix_array = (pos_mix_array + 1) % length_mix_array;
    }
}