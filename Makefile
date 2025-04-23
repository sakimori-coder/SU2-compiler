# ────────────────────────────────────────────
# ディレクトリ定義
SRCDIR := src
INCDIR := include
OBJDIR := obj
BINDIR := bin
TESTDIR := tests


# ────────────────────────────────────────────
# ソース・オブジェクト一覧
SRCS := $(shell find $(SRCDIR) -name '*.cpp')
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

TEST_SRCS  := $(shell find $(TESTDIR) -name '*.cpp')
TEST_OBJS  := $(patsubst $(TESTDIR)/%.cpp,$(OBJDIR)/tests/%.o,$(TEST_SRCS))

DEPS := $(OBJS:.o=.d)

print-%: ; @echo $* = $($*)


# ────────────────────────────────────────────
# コンパイラ設定
CXX = g++
CXXFLAGS = -std=c++20 -O3 -fopenmp -I$(INCDIR) -lquadmath
CXXFLAGS += -MMD -MP

LDLIBS_TEST := -lgtest -lgtest_main -pthread

# PROFILE_FLAGS = -O2 -g -fno-omit-frame-pointer -fno-inline
# CXXFLAGS += $(PROFILE_FLAGS)


# ────────────────────────────────────────────
# obj/ と bin/ の作成
$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/ring:
	mkdir -p $@

$(BINDIR):
	mkdir -p $@


# ────────────────────────────────────────────
# ソース → オブジェクト（同時に .d も生成）
# $< が src/xxx.cpp, $@ が obj/xxx.o
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(OBJDIR)/tests/%.o: $(TESTDIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@



$(BINDIR)/test_runner: $(OBJS) $(TEST_OBJS)
	@mkdir -p $(BINDIR)
	$(CXX) $^ $(CXXFLAGS) $(LDLIBS_TEST) -o $@

.PHONY: test
test: $(BINDIR)/test_runner
	@echo "Running tests…"; ./$(BINDIR)/test_runner



# ────────────────────────────────────────────
# 依存関係ファイルを読み込む
# 存在しなくてもエラーにしない include なので、初回ビルドでも OK
-include $(DEPS)


# ────────────────────────────────────────────
# クリーン
clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean