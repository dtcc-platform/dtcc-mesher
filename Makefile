BUILD_DIR ?= build
UV ?= uv

.PHONY: all install test test-python clean

all:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR)

install:
	$(UV) sync

test-python:
	$(UV) run pytest

test:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	rm -rf $(BUILD_DIR)
