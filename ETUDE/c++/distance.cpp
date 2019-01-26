/* 1. distance
2. max elemant in array
*/
#include <iostream>
#include <iterator>
#include <vector>
 
int main()
{
    std::vector<int> v{3, 1, 4, 2, 0};
    // 1.
    auto distance = std::distance(v.begin(), v.end());
 
    std::cout << distance << '\n';
    // 2.
    int max = std::distance(v.begin(), std::max_element(v.begin(), v.end()));
    std::cout << "max = " << max << std::endl;  
}

