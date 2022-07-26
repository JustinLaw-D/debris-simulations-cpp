#include "Arrays.h"

void array_test() {
    const array<size_t,3> dim = {3,3,3};
    ArrayND<double, 3> arr = ArrayND<double, 3>(0.0, dim);
    arr.print();
    ArrayND<double, 2> * arr2 = new ArrayND<double, 2>(0.0, array<size_t,2>({3,2}));
    arr2->print();
    arr2->at(array<size_t,2>({1,0})) = 10.0;
    arr2->copy_sum(571.0);
    arr2->print();
    ArrayND<bool, 2> arr3 = ArrayND<bool, 2>(true, array<size_t, 2>({6,2}));
    arr3.print();
    ArrayND<bool, 2> * arr4 = new ArrayND<bool, 2>(arr3);
    arr4->print();
    ArrayND<float, 2> * arr5 = new ArrayND<float, 2>(array<size_t, 2>({2,2}));
    for (size_t i = 0; i < 2; i++) {
        for (size_t j = 0; j < 2; j++) {
            arr5->at(array<size_t, 2>({i,j})) = 0.0;
        }
    }
    arr5->print();
}

int main() {
    array_test();
    return 0;
}
