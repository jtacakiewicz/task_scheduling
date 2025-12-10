# Nazwa pliku źródłowego C++
SOURCE = verify.cpp

# Nazwa pliku wykonywalnego
EXECUTABLE = weryfikator

# Kompilator
CXX = g++

# Flagi kompilatora:
# -std=c++11: Użycie standardu C++11 (lub nowszego, np. c++17)
# -Wall: Włączenie wszystkich standardowych ostrzeżeń
# -O2: Poziom optymalizacji
CXXFLAGS = -std=c++17 -Wall -O2

# Cel domyślny: kompilacja programu
all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $@ $^
build: $(EXECUTABLE)

# Cel do czyszczenia skompilowanych plików
clean:
	rm -f $(EXECUTABLE) *.in *.out

# ----------------------------------------------------------------------
# PRZYKŁADY UŻYCIA
# ----------------------------------------------------------------------

# Cel do generowania przykładowej instancji (10 zadań)
generate_ex: $(EXECUTABLE)
	@echo "--- Generowanie przykładowej instancji (10 zadań) do test.in ---"
	./$(EXECUTABLE) generate test.in 5

# Cel do weryfikacji przykładowego wyniku
# UWAGA: Ten cel wymaga, aby plik test.out istniał i zawierał poprawną sekwencję
# oraz oczekiwaną sumę opóźnień.
# Możesz użyć 'generate' jako kroku wstępnego.
verify_ex: $(EXECUTABLE)
	@echo "--- Weryfikacja przykładowego harmonogramu w test.out dla instancji test.in ---"
	@if [ ! -f test.in ]; then \
	    echo "BŁĄD: Plik test.in nie istnieje. Najpierw uruchom 'make generate'."; \
	elif [ ! -f test.out ]; then \
	    echo "OSTRZEŻENIE: Plik test.out nie istnieje. Używam domyślnego pliku z zerowym wynikiem."; \
	    echo "0" > test.out; \
	    echo "1 2 3 4 5" >> test.out; \
	fi
	./$(EXECUTABLE) test.in test.out

.PHONY: all clean generate verify
