# ────────────────────────────────────────────
# ディレクトリ定義
SRCDIR := src
INCDIR := include
TESTDIR := tests
BUILDDIR := build
OBJDIR := $(BUILDDIR)/obj
LIBDIR := $(BUILDDIR)/lib
TEST_BIN := $(BUILDDIR)/test_runner


# ────────────────────────────────────────────
# ソース・オブジェクト一覧
SRCS := $(shell find $(SRCDIR) -name '*.cpp')
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

TEST_SRCS  := $(shell find $(TESTDIR) -name '*.cpp')
TEST_OBJS  := $(patsubst $(TESTDIR)/%.cpp,$(OBJDIR)/tests/%.o,$(TEST_SRCS))

DEPS := $(OBJS:.o=.d) $(TEST_OBJS:.o=.d)

print-%: ; @echo $* = $($*)


# ────────────────────────────────────────────
# コンパイラ設定
CXX = g++
CXXFLAGS = -std=c++20 -O3 -fopenmp -I$(INCDIR) -Iextern/boost -Iextern/eigen
CXXFLAGS += -MMD -MP
CXXFLAGS += -fdiagnostics-show-template-tree
CXXFLAGS += -fno-elide-type
LDFLAGS := -fopenmp
LDLIBS := -pthread -lquadmath

GTEST_INC := -Iextern/googletest/googletest -Iextern/googletest/googletest/include
GTEST_LIB := $(LIBDIR)/libgtest.a
GTEST_MAIN:= $(LIBDIR)/libgtest_main.a

# PROFILE_FLAGS = -O2 -g -fno-omit-frame-pointer -fno-inline
# CXXFLAGS += $(PROFILE_FLAGS)


# ────────────────────────────────────────────
# buildの作成
$(BUILDDIR):
	mkdir $@

$(OBJDIR):
	mkdir -p $@


# ────────────────────────────────────────────
# ソース → オブジェクト（同時に .d も生成）
# $< が src/xxx.cpp, $@ が obj/xxx.o
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(OBJDIR)/tests/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -c $< -o $@



$(TEST_BIN): $(TEST_OBJS) $(GTEST_LIB) $(GTEST_MAIN) $(OBJS)
	@mkdir -p $(dir $@)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

.PHONY: test
test: $(TEST_BIN)
	@echo "Running tests…"; ./$(TEST_BIN)


# ────────────────────────────────────────────
# GoogleTest ライブラリ
GTEST_SRCDIR := extern/googletest/googletest/src

$(GTEST_LIB): $(GTEST_SRCDIR)/gtest-all.cc
	@mkdir -p $(LIBDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -c $< -o $(OBJDIR)/gtest-all.o
	ar rcs $@ $(OBJDIR)/gtest-all.o

$(GTEST_MAIN): $(GTEST_SRCDIR)/gtest_main.cc
	@mkdir -p $(LIBDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -c $< -o $(OBJDIR)/gtest_main.o
	ar rcs $@ $(OBJDIR)/gtest_main.o


# ────────────────────────────────────────────
# 依存関係ファイルを読み込む
# 存在しなくてもエラーにしない include なので、初回ビルドでも OK
-include $(DEPS)


# ────────────────────────────────────────────
# クリーン
clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean