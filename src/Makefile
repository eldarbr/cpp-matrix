CC = g++
CCFL = -Wall -Werror -Wextra --std=c++17
TEST_LIBS = -lgtest

LIB_OBJ_DIR := .obj_lib
TESTS_OBJ_DIR := .obj_tests

all: matrix.a

matrix.a: prepare_l
	rm -f matrix.a
	$(CC) $(CCFL) -c matrix.cpp -o $(LIB_OBJ_DIR)/matrix.o
	ar rcs $@ $(LIB_OBJ_DIR)/matrix.o

test: matrix.a prepare_t
	$(CC) $(CCFL) -c tests/tests.cpp -o $(TESTS_OBJ_DIR)/tests.o
	$(CC) $(CCFL) -o test.out $(LIB_OBJ_DIR)/*.o $(TESTS_OBJ_DIR)/*.o -lm $(TEST_LIBS)
	./test.out

gcov_report: clean prepare_l prepare_t
	$(CC) $(CCFL) --coverage -c matrix.cpp -o $(LIB_OBJ_DIR)/matrix.o
	$(CC) $(CCFL) -c tests/tests.cpp -o $(TESTS_OBJ_DIR)/tests.o
	$(CC) $(CCFL) --coverage $(TESTS_OBJ_DIR)/*.o $(LIB_OBJ_DIR)/*.o -o test_report.out -lm $(TEST_LIBS)
	./test_report.out
	lcov -t test_report -o test.info -c -d . --no-external
	genhtml -o report test.info
	open ./report/index.html

prepare_l:
	mkdir -p $(LIB_OBJ_DIR)

prepare_t:
	mkdir -p $(TESTS_OBJ_DIR)

clean:
	rm -rf .obj* matrix.a test.out
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -f test_report.out
	rm -rf report

rebuild: clean all

format:
	find .. \( -name '*.cpp' -o -name '*.h' \) -print0 | xargs -0 clang-format -n -Werror --verbose --style=google

format-force:
	find .. \( -name '*.cpp' -o -name '*.h' \) -print0 | xargs -0 clang-format -i --verbose --style=google
