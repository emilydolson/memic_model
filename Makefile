# Project-specific settings
PROJECT := memic_model
EMP_DIR := ../Empirical/source

# Flags to use regardless of compiler
CFLAGS_all := -Wall -Wno-unused-function -std=c++17 -I$(EMP_DIR)/

# Native compiler information
CXX_nat := g++
CFLAGS_nat := -O3 -DNDEBUG $(CFLAGS_all)
CFLAGS_nat_debug := -g -DEMP_TRACK_MEM $(CFLAGS_all)

# Emscripten compiler information
CXX_web := emcc
OFLAGS_web_all := -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap']" -s TOTAL_MEMORY=67108864 --js-library $(EMP_DIR)/web/library_emp.js -s EXPORTED_FUNCTIONS="['_main', '_empCppCallback']" -s DISABLE_EXCEPTION_CATCHING=1 -s NO_EXIT_RUNTIME=1 #--embed-file configs
OFLAGS_web := -Oz -DNDEBUG
OFLAGS_web_debug := -g4 -Oz -pedantic -Wno-dollar-in-identifier-extension -s ASSERTIONS=2 -s WASM=0 -s DEMANGLE_SUPPORT=1

CFLAGS_web := $(CFLAGS_all) $(OFLAGS_web) $(OFLAGS_web_all)
CFLAGS_web_debug := $(CFLAGS_all) $(OFLAGS_web_debug) $(OFLAGS_web_all)


default: $(PROJECT)
native: $(PROJECT)
web: $(PROJECT).js
all: $(PROJECT) $(PROJECT).js

debug:	CFLAGS_nat := $(CFLAGS_nat_debug)
debug:	$(PROJECT)

debug-web:	CFLAGS_web := $(CFLAGS_web_debug)
debug-web:	$(PROJECT).js

web-debug:	debug-web

$(PROJECT):	source/native/$(PROJECT).cc
	$(CXX_nat) $(CFLAGS_nat) source/native/$(PROJECT).cc -o $(PROJECT)
	@echo To build the web version use: make web

$(PROJECT).js: source/web/$(PROJECT)-web.cc
	$(CXX_web) $(CFLAGS_web) source/web/$(PROJECT)-web.cc -o web/$(PROJECT).js

test: tests/unit_tests.cc
	$(CXX_nat) $(CFLAGS_nat_debug) tests/unit_tests.cc -o test_debug.out
	./test_debug.out 
	$(CXX_nat) $(CFLAGS_nat) tests/unit_tests.cc -o test_optimized.out
	./test_optimized.out
	make web-debug
	make web

coverage: tests/unit_tests.cc
	cp ../force-cover/force_cover .
	cp ../force-cover/fix_coverage.py .
	rsync -r --exclude .git --exclude web . ../coverage_testing
	tests/convert_for_coverage.sh $(CFLAGS_nat_debug)
	clang++-7 -mllvm -enable-name-compression=false -fprofile-instr-generate -fcoverage-mapping -O0 -fno-inline -fno-elide-constructors $(CFLAGS_nat_debug) ../coverage_testing/tests/unit_tests.cc -o coverage_test.out
	./coverage_test.out
	llvm-profdata merge default.profraw -o default.profdata
	llvm-cov show coverage_test.out -instr-profile=default.profdata > coverage.txt
	python fix_coverage.py coverage.txt
	rm -r ../coverage_testing
	rm force_cover
	rm fix_coverage.py

clean:
	rm -f $(PROJECT) web/$(PROJECT).js web/*.js.map web/*.js.map *~ source/*.o test_debug.out test_optimized.out coverage_test.out coverage.txt default.profdata default.profraw

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
