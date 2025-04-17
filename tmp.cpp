#include <mpi.h>
#include <omp.h>
#include <vector>
#include <tuple>
#include <iostream>

// 関数fのプロトタイプ宣言（例）
int f(int x, int y, int z){
    return x + y + z;
}

// MPI並列化を行うサブルーチン
bool parallel_subroutine(std::vector<int>& X, std::vector<int>& Y, std::vector<int>& Z, int key) {
    MPI_Init(NULL, NULL);  // MPIの初期化

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<std::tuple<int, int, int>> local_results;

    // 各プロセスに配布するXのサブセットを計算
    int chunk_size = X.size() / world_size;
    int start = world_rank * chunk_size;
    int end = (world_rank == world_size - 1) ? X.size() : start + chunk_size;

    // #pragma omp parallel for collapse(2)
    for (int i = start; i < end; ++i) {
#pragma omp parallel for
        for (auto &y : Y) {
#pragma omp parallel for
            for (auto &z : Z) {
                if (f(X[i], y, z) == key) {
                    #pragma omp critical
                    {
                        // 処理
                        std::cout << "Found key at process " << world_rank << ", thread " << omp_get_thread_num() << std::endl;
                        local_results.push_back({X[i], y, z});
                    }
                }
            }
        }
    }

    std::vector<std::tuple<int, int, int>> global_results;
    // Rootプロセスで結果を集約
    if (world_rank == 0) {
        // 各プロセスの結果を受信
        for (int i = 0; i < world_size; ++i) {
            if (i == 0) {
                // 自分のプロセスの結果を追加
                global_results.insert(global_results.end(), local_results.begin(), local_results.end());
            } else {
                // 他のプロセスからの結果を受信
                int recv_size;
                MPI_Recv(&recv_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                std::vector<std::tuple<int, int, int>> recv_results(recv_size);
                MPI_Recv(recv_results.data(), recv_size * sizeof(std::tuple<int, int, int>), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                global_results.insert(global_results.end(), recv_results.begin(), recv_results.end());
            }
        }

        // 結果を表示
        for (const auto &result : global_results) {
            int x, y, z;
            std::tie(x, y, z) = result;
            std::cout << "Result: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }
    } else {
        // 他のプロセスは自分の結果をRootプロセスに送信
        int local_size = local_results.size();
        MPI_Send(&local_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_results.data(), local_size * sizeof(std::tuple<int, int, int>), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();  // MPIの終了

        for (const auto &result : global_results) {
            int x, y, z;
            std::tie(x, y, z) = result;
            std::cout << "Result: (" << x << ", " << y << ", " << z << ")" << std::endl;
        }
    if(!global_results.empty()) return true;
}

int main(int argc, char *argv[]) {
    // X, Y, Zの初期化（例）
    std::vector<int> X = {1, 2, 3, 4, 5};
    std::vector<int> Y = {6, 7, 8, 9, 10};
    std::vector<int> Z = {11, 12, 13, 14, 15};
    int key = 30;

    // MPI並列化を行うサブルーチンの呼び出し
    parallel_subroutine(X, Y, Z, key);

    // 他の処理（MPIを使わない）
    std::cout << "Main function continues without MPI." << std::endl;

    return 0;
}
