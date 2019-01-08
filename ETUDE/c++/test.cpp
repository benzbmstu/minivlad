/*
Содержание:
1) как использовать двумерный вектор, состоящий из пар <int, char>, инициализировать его, 
   а так же сортировать по первому <int> значению. Для этого используется кастомный 
   механизм сортировки с помощью функции compare().
2) как проверить наличие элемента в одномерном векторе, а так же вернуть
   его индекс. 
*/
#include <iostream> 
#include <vector> 
#include <algorithm> 

bool compare(const std::pair<int, int>&i, const std::pair<int, int>&j){ 
    return i.first < j.first; 
} 

int main(){
    /* -1- */ 
    std::vector<std::pair<int, char>> vec; 
    vec.resize(10); // reserve space for 10 elements 
    std::string letters = "abcdefghij"; 
    int randNum; 

    for(int i=0; i<10; i++){ 
        randNum = std::rand()%(10 + 1); // generate random numbers between 0 and 10 
        vec[i].first = randNum; // assign random integer to first element of pair 
        vec[i].second = letters[i]; // assign letter to second element of pair 
    } 
    std::cout << std::fixed;
    std::cout.precision(4);
    for(int i=0; i<10; i++) // print out unsorted array 
        std::cout << vec[i].first << " " << vec[i].second << "\n"; 

    std::cout << "\n"; 
    std::sort(vec.begin(), vec.end(), compare); 

    for(int i=0; i<10; i++) // print out sorted array 
        std::cout << vec[i].first << " " << vec[i].second << "\n"; 

    /* -2- */
    // my tests
    std::vector<int> vecOfNums = { 12, 45, 54, 33, 2, 7, 8, 22, 43, 19 };

    for (int &a : vecOfNums)
        std::cout << a << ' ';
    std::cout << std::endl;

    // Check if element 10 exists in vector
    std::vector<int>::iterator it = std::find(vecOfNums.begin(), vecOfNums.end(), 22);
    if (it != vecOfNums.end())
        std::cout << "Element Found" << std::endl;
    else
        std::cout << "Element Not Found" << std::endl;

    // Get index of element from iterator
    int index = std::distance(vecOfNums.begin(), it);
    std::cout << std::endl << "Его индекс : " << index;

    return 1; 
} 
