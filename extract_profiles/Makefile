CXXFLAGS=-Wall
CC=$(CXX)

SOURCES=extract.c++
TARGET=extract

.PHONY: all
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CC) $< -o $@ -lboost_system -lboost_filesystem

clean:
	rm -rf $(TARGET)
