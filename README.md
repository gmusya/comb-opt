# Dependencies

## OR-Tools
```
git clone https://github.com/google/or-tools
cd or-tools
cmake -S . -B build -DBUILD_DEPS=ON
cmake --build build --config Release --target install -j 16
```
