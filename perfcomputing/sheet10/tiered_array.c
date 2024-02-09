#include "tiered_array.h"

void create_random_element(char* element, int element_size){
    for (int i = 0; i < element_size; i++){
        element[i] = i % 255;
    }
}

void create_random_container(char** container, int element_size){
    for (int i = 0; i < _CONTAINER_SIZE; i++){
        create_random_element(container[i], element_size);
    }
}

tiered_array_t* init_tiered_array(int element_size, int num_elements){
    tiered_array_t* array = malloc(sizeof(tiered_array_t) + element_size * sizeof(char));
    array->num_containers = (num_elements / _CONTAINER_SIZE) + 1;

    // add to num_elements in case there is a num_elements that is 
    // between multiples of _CONTAINER_SIZE
    for (int i = 0; i < num_elements + (_CONTAINER_SIZE -1); i += _CONTAINER_SIZE){
        char** container; 
        for (int j = 0; j < _CONTAINER_SIZE; j++){
            create_random_element(container[j], element_size);
        }
        array->containers[i] = container;
    }
    return array;
}

/* I think either add to current container or create new one,
   just moving the pointers in array->containers */
void tiered_array_insert(tiered_array_t* array, int index, char* value){
    if(array->containers[index/_CONTAINER_SIZE][index % _CONTAINER_SIZE] == NULL){
        array->containers[index / _CONTAINER_SIZE][index % _CONTAINER_SIZE] = value;
    }
    else{
        size_t amount_containers = (sizeof(*array->containers)/sizeof(char***));
        char** new_containers[amount_containers + 1];
        memcpy(new_containers, array->containers, sizeof(*array->containers));
        for (int i = index / _CONTAINER_SIZE; i < (int) amount_containers; i++){
           new_containers[i+1] = new_containers[i];
        }
        new_containers[index / _CONTAINER_SIZE][index % _CONTAINER_SIZE] = value;
        array->containers = new_containers;
    }
} // I am not sure this works, if i could understand this 3 days from now that would be nice

/* remove == set value to Null, see insertion*/
void tiered_array_remove();

void traverse_tiered_array(tiered_array_t* array, int num_elements, int element_size, Op ops[], int num_op){
    
}