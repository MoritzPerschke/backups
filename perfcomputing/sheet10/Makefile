CFLAGS = -std=c11 -Wall -Wextra -g -pedantic -O3 -pthread

SOURCES = benchmark.c array.c linked_list.c unrolled_linked_list.c tiered_array.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = benchmark

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

array.o: array.h
linked_list.o: linked_list.h

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJECTS)