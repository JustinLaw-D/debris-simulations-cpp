// 1 and N-d, constant size arrays

#include <cstring>
#include <cmath>
#include <iostream>
#include <array>
#pragma once

using namespace std;

// classes for constant-size N-d arrays (as well as the 1-d subclass)

template <typename T, size_t N>
class ArrayND
{
    protected:
        T * arr; // pointer to the array
        array<size_t, N> dim; // list of dimentions
        size_t tot_size; // total number of elements
        ArrayND(T * arr, const array<size_t, N> dim, const size_t tot_size); // explicit constructor
        void print_dim(array<size_t, N> partial_loc, size_t num_set, size_t curr_true_loc, size_t curr_skip_mul);

    public:
        ArrayND(); // basic constructor, creates size zero array
        ArrayND(const array<size_t, N> &dim); // another basic constructor 
        ArrayND(const ArrayND &arr); // copy constructor

        static ArrayND zeroes_int(const array<size_t, N> &dim) {
            /*
            initializes an array of all zeros, of type int

            Input(s):
            dim : list of dimentions of the array

            Output(s):
            arr_out : ArrayND object full of zeros
            */
            if (N == 0) { // in this case, return zero array
                return ArrayND<int, N>();
            } else {
                size_t tot_size = 1; // calculate total number of elements in the array
                for (size_t i = 0; i < N; i++) {
                    tot_size *= dim[i];
                }
                int * arr = new int[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0;
                }
                return ArrayND<int, N>(arr, dim_loc, tot_size);
            }
        }

        static ArrayND * zeroes_int_dyn(const array<size_t, N> &dim) {
            /*
            initializes a dynamically allocated array of all zeros, of type int

            Input(s):
            dim : list of dimentions of the array

            Output(s):
            arr_out : ArrayND object full of zeros
            */
            if (N == 0) { // in this case, return zero array
                return new ArrayND<int, N>();
            } else {
                size_t tot_size = 1; // calculate total number of elements in the array
                for (size_t i = 0; i < N; i++) {
                    tot_size *= dim[i];
                }
                int * arr = new int[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0;
                }
                return new ArrayND<int, N>(arr, dim_loc, tot_size);
            }
        }

        static ArrayND zeroes_double(const array<size_t, N> &dim) {
            /*
            initializes a dynamically allocated array of all zeros, of type int

            Input(s):
            dim : list of dimentions of the array

            Output(s):
            arr_out : ArrayND object full of zeros
            */
            if (N == 0) { // in this case, return zero array
                return ArrayND<double, N>();
            } else {
                size_t tot_size = 1; // calculate total number of elements in the array
                for (size_t i = 0; i < N; i++) {
                    tot_size *= dim[i];
                }
                double * arr = new double[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0.0;
                }
                return ArrayND<double, N>(arr, dim_loc, tot_size);
            }
        }

        static ArrayND * zeroes_double_dyn(const array<size_t, N> &dim) {
            /*
            initializes a dynamically allocated array of all zeros, of type double

            Input(s):
            dim : list of dimentions of the array
            num_dim : number of dimentions of the array

            Output(s):
            arr_out : ArrayND object full of zeros
            */
            if (N == 0) { // in this case, return zero array
                return new ArrayND<double, N>();
            } else {
                size_t tot_size = 1; // calculate total number of elements in the array
                for (size_t i = 0; i < N; i++) {
                    tot_size *= dim[i];
                }
                double * arr = new double[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0.0;
                }
                return new ArrayND<double, N>(arr, dim_loc, tot_size);
            }
        }

        static ArrayND fill(const T val, const array<size_t, N> &dim); // creates array with copies of given value
        static ArrayND * fill_dyn(const T val, const array<size_t, N> &dim); // same thing, but returns a pointer
        inline bool is_zero(); // check if this is a size-zero array
        inline array<size_t,N> get_dim() const; // gets dimensions
        inline size_t get_tot_size() const; // get total number of elements
        inline T * get_arr() const; // gets the underlying array
        inline T & at(const array<size_t, N> &loc); // get reference to a specific element
        void copy_sum(const T num); // basic arithmatic with constant values, modifying current array
        void copy_sub(const T num);
        void copy_mul(const T num);
        void copy_div(const T num);
        void copy_pow(const T num);
        void print(); // displays the array
        ~ArrayND(); // destructor
};

// class for 1D arrays, with some speed shortcuts

template<typename T>
class Array1D : public ArrayND<T,1> {

    private:
        Array1D(T * arr, const array<size_t, 1> dim, const size_t tot_size); // explicit constructor

    public:
        Array1D(); // basic constructor, creates size zero array
        Array1D(const size_t len); // another basic constructor 
        Array1D(const Array1D &arr); // copy constructor

        static Array1D zeroes_int(const size_t len) {
            /*
            initializes an array of all zeros, of type int

            Input(s):
            len : length of the array

            Output(s):
            arr_out : Array1D object full of zeros
            */
            if (len == 0) { // in this case, return zero array
                return Array1D<int>();
            } else {
                size_t tot_size = len; // calculate total number of elements in the array
                int * arr = new int[tot_size]; array<size_t,1> dim_loc = array<size_t,1>({len});
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0;
                }
                return Array1D<int>(arr, dim_loc, tot_size);
            }
        }

        static Array1D * zeroes_int_dyn(const size_t len) {
            /*
            initializes a dynamically allocated array of all zeros, of type int

            Input(s):
            len : length of the array

            Output(s):
            arr_out : Array1D object full of zeros
            */
            if (len == 0) { // in this case, return zero array
                return new Array1D<int>();
            } else {
                size_t tot_size = len; // calculate total number of elements in the array
                int * arr = new int[tot_size]; array<size_t,1> dim_loc = array<size_t,1>({len});
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0;
                }
                return new Array1D<int>(arr, dim_loc, tot_size);
            }
        }

        static Array1D zeroes_double(const size_t len) {
            /*
            initializes an array of all zeros, of type double

            Input(s):
            len : length of the array

            Output(s):
            arr_out : Array1D object full of zeros
            */
            if (len == 0) { // in this case, return zero array
                return Array1D<double>();
            } else {
                size_t tot_size = len; // calculate total number of elements in the array
                double * arr = new double[tot_size]; array<size_t,1> dim_loc = array<size_t,1>({len});
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0.0;
                }
                return Array1D<double>(arr, dim_loc, tot_size);
            }
        }

        static Array1D * zeroes_double_dyn(const size_t len) {
            /*
            initializes a dyanamically allocated array of all zeros, of type double

            Input(s):
            len : length of the array

            Output(s):
            arr_out : Array1D object full of zeros
            */
            if (len == 0) { // in this case, return zero array
                return new Array1D<double>();
            } else {
                size_t tot_size = len; // calculate total number of elements in the array
                double * arr = new double[tot_size]; array<size_t,1> dim_loc = array<size_t,1>({len});
                // initialize to zeros
                for (size_t i = 0; i < tot_size; i++) {
                    arr[i] = 0.0;
                }
                return new Array1D<double>(arr, dim_loc, tot_size);
            }
        }

        static Array1D fill(const T val, const size_t len); // creates array with copies of given value
        static Array1D * fill_dyn(const T val, const size_t len); // same thing, but returns a pointer
        inline T & at(const size_t loc); // get reference to a specific element
};

template <typename T, size_t N>
ArrayND<T, N>::ArrayND() {
    /*
    creates size-zero ND array

    Input(s): None

    Output(s):
    arr_out : ArrayND with unititialized array/dimentions and size zero
    */
    this->tot_size = 0;
}

template <typename T, size_t N>
ArrayND<T, N>::ArrayND(const array<size_t, N> &dim) {
    /*
    most basic initializer for a n-d array

    Input(s):
    dim : list of dimentions of the array

    Output(s):
    arr_out : ArrayND object with specified size
    */
    size_t tot_size = 1; // calculate total number of elements in the array
    for (size_t i = 0; i < N; i++) {
        tot_size *= dim[i];
    }
    this->arr = new T[tot_size]; this->dim = array<size_t, N>(dim); this->tot_size = tot_size;
}

template <typename T, size_t N>
ArrayND<T, N>::ArrayND(T * arr, const array<size_t, N> dim, const size_t tot_size) {
    /*
    most basic initializer for a n-d array

    Input(s):
    arr : nd array of values
    dim : list of dimentions of the array
    tot_size : the total number of elements in the array

    Output(s):
    arr_out : ArrayND object with specified parameters

    Note(s): expects to take ownership of all pointers given
    */
    this->arr = arr; this->dim = array<size_t, N>(dim); this->tot_size = tot_size;
}

template <typename T, size_t N>
ArrayND<T, N>::ArrayND(const ArrayND<T, N> &arr) {
    /*
    copy constructor

    Input(s):
    arr : array to copy

    Output(s):
    arr_out : deep copy of given array
    */ 
    size_t tot_size = arr.get_tot_size();
    if (tot_size == 0) {this->tot_size = 0;} // size zero array case
    else {
        this->arr = new T[tot_size]; // create new array
        memcpy(this->arr, arr.get_arr(), sizeof(T)*tot_size); // copy over contents of old array
        this->dim = array<size_t, N>(arr.get_dim()); this->tot_size = tot_size;
    }
}

template <typename T, size_t N>
ArrayND<T,N> ArrayND<T,N>::fill(const T val, const array<size_t, N> &dim) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    dim : dimentions of the array

    Output(s):
    arr_out : ArrayND object full of specified value
    */
    if (N == 0) { // in this case, return zero array
        return ArrayND<T, N>();
    } else {
        size_t tot_size = 1; // calculate total number of elements in the array
        for (size_t i = 0; i < N; i++) {
            tot_size *= dim[i];
        }
        T * arr = new T[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
        // initialize to zeros
        for (size_t i = 0; i < tot_size; i++) {
            arr[i] = val;
        }
        return ArrayND<T, N>(arr, dim_loc, tot_size);
    }
}

template <typename T, size_t N>
ArrayND<T, N> * ArrayND<T, N>::fill_dyn(const T val, const array<size_t, N> &dim) {
    /*
    initializes a pointer to an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    dim : dimentions of the array

    Output(s):
    arr_out : pointer to Array2D object full of specified value
    */
    if (N == 0) { // in this case, return zero array
        return new ArrayND<T, N>();
    } else {
        size_t tot_size = 1; // calculate total number of elements in the array
        for (size_t i = 0; i < N; i++) {
            tot_size *= dim[i];
        }
        T * arr = new T[tot_size]; array<size_t,N> dim_loc = array<size_t,N>(dim);
        // initialize to zeros
        for (size_t i = 0; i < tot_size; i++) {
            arr[i] = val;
        }
        return new ArrayND<T, N>(arr, dim_loc, tot_size);
    }
}

template <typename T, size_t N>
array<size_t, N> ArrayND<T, N>::get_dim() const {
    /*
    returns the dimensions of the array
    */
    return this->dim;
}

template <typename T, size_t N>
size_t ArrayND<T, N>::get_tot_size() const {
    /*
    returns the underlying array pointer, should not be used
    */
    return this->tot_size;
}

template <typename T, size_t N>
T * ArrayND<T, N>::get_arr() const {
    /*
    returns the underlying array pointer, should not be used
    */
    return this->arr;
}

template <typename T, size_t N>
bool ArrayND<T, N>::is_zero() {
    /*
    returns true on a size zero array
    */
    return (this->tot_size == 0);
}

template <typename T, size_t N>
T & ArrayND<T, N>::at(const array<size_t,N> &loc) {
    /*
    gets reference to given location in the array

    Input(s):
    loc : location in the array

    Output(s):
    val : reference to the array at that location, performs no size check
    */
    size_t real_loc = 0;
    size_t skip_mul = 1;
    for (size_t i = 0; i < N; i++) {
        real_loc += loc[i]*skip_mul;
        skip_mul *= (this->dim)[i];
    }
    return (this->arr)[real_loc];
}

template <typename T, size_t N>
void ArrayND<T, N>::copy_sum(const T num) {
    /*
    basic addition of a constant, overwriting current matrix
    
    Input(s):
    num : value to be added

    Output(s): None
    */
    for (size_t i = 0; i < this->tot_size; i++) {
        (this->arr)[i] += num;
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::copy_sub(const T num) {
    /*
    basic subtraction of a constant, overwriting current matrix
    
    Input(s):
    num : value to be subtracted

    Output(s): None
    */
    for (size_t i = 0; i < this->tot_size; i++) {
        (this->arr)[i] -= num;
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::copy_mul(const T num) {
    /*
    basic multiplication by a constant, overwriting current matrix
    
    Input(s):
    num : value to be multiplied by

    Output(s): None
    */
    for (size_t i = 0; i < this->tot_size; i++) {
        (this->arr)[i] *= num;
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::copy_div(const T num) {
    /*
    basic division by a constant, overwriting current matrix
    
    Input(s):
    num : value to be divided by

    Output(s): None
    */
    for (size_t i = 0; i < this->tot_size; i++) {
        (this->arr)[i] /= num;
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::copy_pow(const T num) {
    /*
    basic exponentiation by a constant, overwriting current matrix
    
    Input(s):
    num : value to be exponentiated by

    Output(s): None
    */
    for (size_t i = 0; i < this->tot_size; i++) {
        (this->arr)[i] = pow((this->arr)[i], num);
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::print() {
    // prints out array
    if (N == 0) {
        cout << "Zero-Sized Array" << endl;
    } else {
        array<size_t, N> partial_loc = array<size_t, N>(this->dim);
        cout << "[";
        for (size_t i = 0; i < (this->dim)[0]; i++) {
            partial_loc[0] = i;
            this->print_dim(partial_loc, 1, i, (this->dim)[0]);
        }
        cout << "]" << endl;
    }
}

template <typename T, size_t N>
void ArrayND<T, N>::print_dim(array<size_t, N> partial_loc, size_t num_set, size_t curr_true_loc, size_t curr_skip_mul) {
    /*
    local array printing function, printing a specific dimention/set of them

    Input(s):
    partial_loc : array with specified dimensions to print
    num_set : number of specified dimensions
    curr_true_loc : current running tally of true location
    curr_skip_mul : current skip multiplier

    Output(s) : None
    */
    if (num_set == N) { // print a value
        if (partial_loc[N-1] == (this->dim[N-1]) - 1) { // check if it's the end of the row
            cout << (this->arr)[curr_true_loc];
        } else {
            cout << (this->arr)[curr_true_loc] << ",";
        }
    } else { // continue the recursion
        cout << "[";
        for (size_t i = 0; i < (this->dim)[num_set]; i++) {
            partial_loc[num_set] = i;
            this->print_dim(partial_loc, num_set + 1, curr_true_loc + i*curr_skip_mul, curr_skip_mul*(this->dim)[num_set]);
        }
        if (partial_loc[num_set-1] == (this->dim[num_set-1]) - 1) {cout << "]";} // check if it's the end of the row
        else {cout << "]" << endl;}
    }
}

template <typename T, size_t N>
ArrayND<T, N>::~ArrayND() {
    // destructor for ND array
    delete [] this->arr;
}

template <typename T>
Array1D<T>::Array1D() {
    /*
    creates size-zero 1D array

    Input(s): None

    Output(s):
    arr_out : Array1D with unititialized array/dimentions and size zero
    */
    this->tot_size = 0;
}

template <typename T>
Array1D<T>::Array1D(const size_t len) {
    /*
    most basic initializer for a 1-d array

    Input(s):
    len : length of the array

    Output(s):
    arr_out : Array1D object with specified size
    */
    this->arr = new T[len]; this->dim = array<size_t, 1>({len}); this->tot_size = len;
}

template <typename T>
Array1D<T>::Array1D(T * arr, const array<size_t, 1> dim, const size_t tot_size) {
    /*
    most basic initializer for a n-d array

    Input(s):
    arr : nd array of values
    dim : list of dimentions of the array
    tot_size : the total number of elements in the array

    Output(s):
    arr_out : Array1D object with specified parameters

    Note(s): expects to take ownership of all pointers given
    */
    this->arr = arr; this->dim = array<size_t, 1>(dim); this->tot_size = tot_size;
}

template <typename T>
Array1D<T>::Array1D(const Array1D<T> &arr) {
    /*
    copy constructor

    Input(s):
    arr : array to copy

    Output(s):
    arr_out : deep copy of given array
    */ 
    size_t tot_size = arr.get_tot_size();
    if (tot_size == 0) {this->tot_size = 0;} // size zero array case
    else {
        this->arr = new T[tot_size]; // create new array
        memcpy(this->arr, arr.get_arr(), sizeof(T)*tot_size); // copy over contents of old array
        this->dim = array<size_t, 1>(arr.get_dim()); this->tot_size = tot_size;
    }
}

template <typename T>
Array1D<T> Array1D<T>::fill(const T val, const size_t len) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    len : length of the array

    Output(s):
    arr_out : Array1D object full of specified value
    */
    if (len == 0) { // in this case, return zero array
        return Array1D<T>();
    } else {
        T * arr = new T[len]; array<size_t,1> dim_loc = array<size_t,1>({len});
        // initialize to zeros
        for (size_t i = 0; i < len; i++) {
            arr[i] = val;
        }
        return Array1D<T>(arr, dim_loc, len);
    }
}

template <typename T>
Array1D<T> * Array1D<T>::fill_dyn(const T val, const size_t len) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    len : length of the array

    Output(s):
    arr_out : Array1D object full of specified value
    */
    if (len == 0) { // in this case, return zero array
        return new Array1D<T>();
    } else {
        T * arr = new T[len]; array<size_t,1> dim_loc = array<size_t,1>({len});
        // initialize to zeros
        for (size_t i = 0; i < len; i++) {
            arr[i] = val;
        }
        return new Array1D<T>(arr, dim_loc, len);
    }
}

template <typename T>
T & Array1D<T>::at(const size_t loc) {
    /*
    gets reference to the given element of the array

    Input(s):
    loc : location in the array

    Output(s):
    ref : reference to desired element

    Note(s): performs no checks on the index asked for
    */
    return (this->arr)[loc];
}