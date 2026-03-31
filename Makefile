BUILD_DIR ?= build

.PHONY: all test clean

all:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR)

test:
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR)
	ctest --test-dir $(BUILD_DIR) --output-on-failure

clean:
	rm -rf $(BUILD_DIR)
