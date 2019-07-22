set(GCC_CXX_FLAGS_DEBUG " -pedantic -ftemplate-depth=1024 -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unused")

set(CLANG_CXX_FLAGS_DEBUG " -Weverything -Wno-c++98-compat -Wno-documentation -Wno-unknown-pragmas -Wno-global-constructors -Wno-exit-time-destructors -Wno-padded -Wno-c++98-compat-pedantic -Wno-disabled-macro-expansion")

set(MSVC_CXX_FLAGS_DEBUG " /W4 /FS /MP ")
set(MSVC_CXX_FLAGS_RELEASE " /MP ")
