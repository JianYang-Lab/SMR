.DEFAULT_GOAL := help
.PHONY: help build fmt fmt-check clean distclean

SOURCES := $(wildcard src/*.cpp src/*.hpp)

help: ## Show this help message
	@echo "SMR Makefile targets:"
	@echo ""
	@perl -ne 'printf "  make %-14s %s\n", $$1, $$2 if /^([a-zA-Z_-]+):.*?##\s+(.*)/' $(MAKEFILE_LIST)

build: ## Build smr (runs CMake configure if needed)
	./scripts/local/build.sh -g

fmt: ## Format sources (clang-format + #pragma omp re-indent)
	clang-format -i $(SOURCES)
	perl scripts/indent_omp_pragma.pl $(SOURCES)

fmt-check: ## Check formatting without modifying files
	@fail=0; \
	for f in $(SOURCES); do \
		tmp=$$(mktemp); \
		clang-format "$$f" > "$$tmp" && perl scripts/indent_omp_pragma.pl "$$tmp" >/dev/null; \
		if ! cmp -s "$$f" "$$tmp"; then echo "Would reformat: $$f"; fail=1; fi; \
		rm -f "$$tmp"; \
	done; \
	if [ $$fail -ne 0 ]; then echo "Run 'make fmt' to fix."; exit 1; fi; \
	echo "All sources are properly formatted."

clean: ## Clean compiled objects (keeps CMake cache and fetched deps)
	@if [ -d build/RelWithDebInfo ]; then cmake --build build/RelWithDebInfo --target clean; \
	else echo "Nothing to clean"; fi

distclean: ## Remove all build directories
	rm -rf build
