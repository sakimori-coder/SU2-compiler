#include <iostream>
#include <unordered_map>

// 独自のキー型
namespace SU2_Compiler {
    struct ZRoot2 {
        int x;
        int y;

        // 等値比較のためにoperator==を定義
        bool operator==(const ZRoot2& other) const {
            return x == other.x && y == other.y;
        }
    };
}

// MyKeyのハッシュ関数を定義
namespace std {
    template <>
    struct hash<SU2_Compiler::ZRoot2> {
        std::size_t operator()(const SU2_Compiler::ZRoot2& key) const {
            // ハッシュ関数を定義（簡単な例）
            return std::hash<int>()(key.x) ^ (std::hash<int>()(key.y) << 1);
        }
    };
}

int main() {
    // std::unordered_mapにZRoot2をキーとして使用
    std::unordered_map<SU2_Compiler::ZRoot2, std::string> myMap;

    // 要素を挿入
    myMap[{1, 2}] = "Point A";
    myMap[{3, 4}] = "Point B";

    // 要素を取得して表示
    std::cout << "({1, 2}): " << myMap[{1, 2}] << std::endl;
    std::cout << "({3, 4}): " << myMap[{3, 4}] << std::endl;

    return 0;
}
