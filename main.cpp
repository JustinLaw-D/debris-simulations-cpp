// this is just here for basic testing

#include "Arrays.h"
#include <iostream>
using namespace std;

int main() {

    Array2D<double> arr = Array2D<double>::zeroes_double(3, 3);
    arr.print();
    cout << endl;
    arr.set(10.0,1,2);
    arr.print();
    cout << endl;
    Array2D<double> arr2 = Array2D<double>(arr);
    arr2.set(50.0,0,0);
    arr2.print();
    cout << endl;
    arr.print();
    cout << endl;
    arr2.get_range(0, 2, 0, 2).print();
    cout << endl;
    Array2D<double> sum = arr + arr2;
    sum.print();
    cout << endl;
    arr.copy_sum(arr2);
    arr.print();
    cout << endl;
    Array2D<int> vec = Array2D<int>::zeroes_int(1,5);
    vec.print();
    cout << endl;
    vec.swap_dim();
    vec.print();
    return 0;
}