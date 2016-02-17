CXXFLAGS=-Wall
CC=$(CXX)

SOURCES=extract.c++
TARGET=extract

.PHONY: all
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $< -o $@

clean:
	rm -rf $(TARGET)
