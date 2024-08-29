#pragma once

#include <vector>
#include <functional>
#include <mutex>


namespace SU2_Compiler{
    template <typename KeyType, typename ValueType>
    class HashTable
    {
    private:
        /* data */
        std::vector<ValueType> *table;
        std::size_t size;
        std::hash<KeyType> h;
        std::mutex mtx;
    public:
        HashTable(std::size_t size){
            this->size = size;
            this->table = new std::vector<ValueType>[size];
        }
        ~HashTable(){
            delete [] table;
        }

        inline void insert(const KeyType& key, const ValueType& value){
            std::lock_guard<std::mutex> lock(mtx);
            this->table[h(key) % size].push_back(value);
        }

        inline std::vector<ValueType>& find(const KeyType& key) const {
            return table[h(key) % size];
        }
    };
}