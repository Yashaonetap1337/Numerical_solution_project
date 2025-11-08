CC = g++
CFLAGS = -c -std=c++17 -O2 -Wall -Wextra -pedantic
LDFLAGS =
OBJDIR = build
EXECUTABLE = solver

# Находим все .cpp файлы в текущей директории
SOURCES = $(wildcard *.cpp)
# Преобразуем пути *.cpp → build/*.o
OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SOURCES))

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

# Создание директории build при необходимости
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Правило для сборки .o файлов в build/
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJDIR) $(EXECUTABLE)

distclean: clean
	rm -f numerical_solution.csv analytical_solution.csv

rebuild: clean all

print:
	@echo "Sources: $(SOURCES)"
	@echo "Objects: $(OBJECTS)"

.PHONY: all clean distclean rebuild print