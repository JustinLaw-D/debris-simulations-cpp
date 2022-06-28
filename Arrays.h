// 1 and 2-d, constant size arrays

#include <cstring>
#include <cmath>
#include <iostream>
#pragma once

using namespace std;

// classes for constant-size 1/2-d arrays

template <typename T>
class Array2D
{
    private:
        T * arr; // pointer to the array
        size_t rows; // number of rows and columns
        size_t cols;

    public:
        Array2D(); // basic constructor, creates size zero array
        Array2D(T * arr, const size_t rows, const size_t cols); // explicit constructor
        Array2D(const Array2D &arr); // copy constructor
        static Array2D zeroes_int(const size_t rows, const size_t cols); // creates array of all zeros
        static Array2D zeroes_double(const size_t rows, const size_t cols); // creates array of all zeros
        static Array2D fill(const T val, const size_t rows, const size_t cols); // creates array with copies of given value
        inline bool is_zero(); // check if this is a size-zero array
        inline size_t get_rows() const; // get the number of rows/columns
        inline size_t get_cols() const;
        inline T * get_arr() const; // gets the underlying array, shouldn't generally be used
        inline T get(const size_t row, const size_t col) const; // basic getter
        inline bool set(const T val, const size_t row, const size_t col); // basic setter
        Array2D get_range(const size_t row_min, const size_t row_max, 
                          const size_t col_min, const size_t col_max);
        Array2D operator + (Array2D const &arr); // basic addition
        bool copy_sum(Array2D const &arr); // basic addition, writing to current object
        Array2D operator - (Array2D const &arr); // basic subtraction
        bool copy_sub(Array2D const &arr); // basic subtraction, writing to current object
        Array2D operator * (Array2D const &arr); // basic multiplication
        bool copy_mul(Array2D const &arr); // basic multiplication, writing to current object
        Array2D operator / (Array2D const &arr); // basic division
        bool copy_div(Array2D const &arr); // basic division, writing to current object
        Array2D pow(const T num); // raises array to a constant power
        void copy_pow(const T num); // raises array to a constant power, overwriting current object
        bool swap_dim(); // used for switching from row to column vector
        void print(); // displays the array
        ~Array2D(); // destructor
};

template <typename T>
Array2D<T>::Array2D() {
    /*
    creates size-zero 2D array

    Input(s): None

    Output(s):
    arr_out : Array2D with unititialized array and size zero
    */
    this->rows = 0; this->cols = 0;
}

template <typename T>
Array2D<T>::Array2D(T * arr, const size_t rows, const size_t cols) {
    /*
    most basic initializer for a 2-d array

    Input(s):
    arr : 2d array of values
    rows : number of rows in the array
    cols : number of columns in the array

    Output(s):
    arr_out : Array2D object with specified parameters
    */
    this->arr = arr; this->rows = rows; this->cols = cols;
}

template <typename T>
Array2D<T>::Array2D(const Array2D<T> &arr) {
    /*
    copy constructor

    Input(s):
    arr : array to copy

    Output(s):
    arr_out : deep copy of given array
    */
    size_t rows = arr.get_rows(); // get the array size
    size_t cols = arr.get_cols(); 
    if (rows*cols== 0) {this->rows = 0; this->cols = 0;}
    else {
        this->arr = new T[rows*cols]; // create new array
        memcpy(this->arr, arr.get_arr(), sizeof(T)*rows*cols); // copy over contents of old array
        this->rows = rows; this->cols = cols;
    }
}

template <>
Array2D<int> Array2D<int>::zeroes_int(const size_t rows, const size_t cols) {
    /*
    initializes an array of all zeros, of type int

    Input(s):
    rows : number of rows in the array
    cols : number of columns in the array

    Output(s):
    arr_out : Array2D object full of zeros
    */
    int * arr = new int[rows*cols];
    for (size_t i = 0; i < rows*cols; i++) {arr[i] = 0;}
    return Array2D<int>(arr, rows, cols);
}

template <>
Array2D<double> Array2D<double>::zeroes_double(const size_t rows, const size_t cols) {
    /*
    initializes an array of all zeros, of type double

    Input(s):
    rows : number of rows in the array
    cols : number of columns in the array

    Output(s):
    arr_out : Array2D object full of zeros
    */
    double * arr = new double[rows*cols];
    for (size_t i = 0; i < rows*cols; i++) {arr[i] = 0.0;}
    return Array2D<double>(arr, rows, cols);
}

template <typename T>
Array2D<T> Array2D<T>::fill(const T val, const size_t rows, const size_t cols) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    rows : number of rows in the array
    cols : number of columns in the array

    Output(s):
    arr_out : Array2D object full of specified value
    */
    T * arr = new T[rows*cols];
    for (size_t i = 0; i < rows*cols; i++) {arr[i] = val;}
    return Array2D<T>(arr, rows, cols);
}

template <typename T>
size_t Array2D<T>::get_rows() const {
    /*
    returns the number of rows in the array
    */
    return this->rows;
}

template <typename T>
size_t Array2D<T>::get_cols() const {
    /*
    returns the number of columns in the array
    */
    return this->cols;
}

template <typename T>
T * Array2D<T>::get_arr() const {
    /*
    returns the underlying array pointer, should not be used
    */
    return this->arr;
}

template <typename T>
bool Array2D<T>::is_zero() {
    /*
    returns true on a size zero array
    */
    if ((this->rows)*(this->cols) == 0) {return true;}
    else {return false;}
}

template <typename T>
T Array2D<T>::get(const size_t row, const size_t col) const {
    /*
    basic getter method

    Input(s):
    row : row of the element to get
    col : column of the element to get

    Output(s):
    val : value of the array at that location, or at (0,0) if the
          location is out of bounds
    */
    if ((row >= this->rows) || (col >= this->cols)) {return (this->arr)[0];}
    else {return (this->arr)[row*(this->cols) + col];}
}

template <typename T>
bool Array2D<T>::set(const T val, const size_t row, const size_t col) {
    /*
    basic setter method

    Input(s):
    val : value to set the element to
    row : row of the element to set
    col : column of the element to set

    Output(s):
    suc : true if given indices are in bounds, false otherwise
    */
    if ((row >= this->rows) || (col >= this->cols)) {return false;}
    else {(this->arr)[row*(this->cols) + col] = val; return true;}
}

template <typename T>
Array2D<T> Array2D<T>::get_range(const size_t row_min, const size_t row_max, 
                                 const size_t col_min, const size_t col_max) {
    /*
    gets parts of array within given range

    Input(s):
    row_min : minimum row to get
    row_max : maximum row to get
    col_min : minimum column to get
    col_max : maximum column to get

    Output(s):
    arr_out : array containing the given range, or size zero array if range is invalid

    Note(s): range is in the form [min, max)
    */
    if ((row_min >= this->rows) || (row_max >= this->rows)) {return Array2D<T>();} // check ranges
    if ((col_min >= this->cols) || (row_max >= this->cols)) {return Array2D<T>();}
    size_t num_rows = row_max - row_min; // get size of return array, and create it
    size_t num_cols = col_max - col_min;
    T * arr_new = new T[num_rows*num_cols];
    for (size_t i = 0; i < num_rows; i++) { // iterate through and fill the array
        for (size_t j = 0; j < num_cols; j++) {
            arr_new[i*num_cols + j] = (this->arr)[(row_min + i)*(this->cols) + col_min + j];
        }
    }
    return Array2D<T>(arr_new, num_rows, num_cols);
}

template <typename T>
Array2D<T> Array2D<T>::operator + (Array2D<T> const &arr) {
    /*
    basic addition between matrices
    
    Input(s):
    arr : matrix being added to the current matrix

    Output(s):
    sum : sum of the two matrices, or zero-size matrix if they are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return Array2D<T>();} // handle case of different sizes
    else {
        Array2D<T> sum = Array2D<T>(arr); // get copy of given array
        T * sum_arr = sum.get_arr(); // get the underlying array
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                sum_arr[i*(this->cols) + j] += (this->arr)[i*(this->cols) + j];
            }
        }
        return sum;
    }
}

template <typename T>
bool Array2D<T>::copy_sum(Array2D<T> const &arr) {
    /*
    basic addition between matrices, overwriting the current matrix
    
    Input(s):
    arr : matrix being added to the current matrix

    Output(s):
    suc : true if matrices are the same size, false otherwise

    Note(s): does nothing if matrices are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return false;} // handle case of different sizes
    else {
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                (this->arr)[i*(this->cols) + j] += arr.get(i, j);
            }
        }
        return true;
    }
}

template <typename T>
Array2D<T> Array2D<T>::operator - (Array2D<T> const &arr) {
    /*
    basic subtraction between matrices
    
    Input(s):
    arr : matrix being subtracted from the current matrix

    Output(s):
    diff : difference of the two matrices, or zero-size matrix if they are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return Array2D<T>();} // handle case of different sizes
    else {
        Array2D<T> diff = Array2D<T>(this); // get copy of given array
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                diff.set(diff.get(i,j) - arr.get(i,j), i, j);
            }
        }
        return diff;
    }
}

template <typename T>
bool Array2D<T>::copy_sub(Array2D<T> const &arr) {
    /*
    basic subtraction between matrices, overwriting the current matrix
    
    Input(s):
    arr : matrix being subtracted from the current matrix

    Output(s):
    suc : true if matrices are the same size, false otherwise

    Note(s): does nothing if matrices are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return false;} // handle case of different sizes
    else {
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                (this->arr)[i*(this->cols) + j] -= arr.get(i, j);
            }
        }
        return true;
    }
}

template <typename T>
Array2D<T> Array2D<T>::operator * (Array2D<T> const &arr) {
    /*
    basic multiplication between matrices
    
    Input(s):
    arr : matrix being multiplied with the current matrix

    Output(s):
    prod : product of the two matrices, or zero-size matrix if they are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return Array2D<T>();} // handle case of different sizes
    else {
        Array2D<T> prod = Array2D<T>(arr); // get copy of given array
        T * prod_arr = prod.get_arr(); // get the underlying array
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                prod_arr[i*(this->cols) + j] *= (this->arr)[i*(this->cols) + j];
            }
        }
        return prod;
    }
}

template <typename T>
bool Array2D<T>::copy_mul(Array2D<T> const &arr) {
    /*
    basic multiplication between matrices, overwriting the current matrix
    
    Input(s):
    arr : matrix being multiplied with the current matrix

    Output(s):
    suc : true if matrices are the same size, false otherwise

    Note(s): does nothing if matrices are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return false;} // handle case of different sizes
    else {
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                (this->arr)[i*(this->cols) + j] *= arr.get(i, j);
            }
        }
        return true;
    }
}

template <typename T>
Array2D<T> Array2D<T>::operator / (Array2D<T> const &arr) {
    /*
    basic division between matrices
    
    Input(s):
    arr : matrix dividing the current matrix

    Output(s):
    prod : quotient of the two matrices, or zero-size matrix if they are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return Array2D<T>();} // handle case of different sizes
    else {
        Array2D<T> prod = Array2D<T>(this); // get copy of given array
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                prod.set(prod.get(i,j)/arr.get(i,j), i, j)
            }
        }
        return prod;
    }
}

template <typename T>
bool Array2D<T>::copy_div(Array2D<T> const &arr) {
    /*
    basic division between matrices, overwriting the current matrix
    
    Input(s):
    arr : matrix being divided into the current matrix

    Output(s):
    suc : true if matrices are the same size, false otherwise

    Note(s): does nothing if matrices are different sizes
    */
    if ((arr.get_rows() != this->rows) || (arr.get_cols() != this->cols)) {return false;} // handle case of different sizes
    else {
        for (size_t i = 0; i < this->rows; i++) {
            for (size_t j = 0; j < this->cols; j++) {
                (this->arr)[i*(this->cols) + j] /= arr.get(i, j);
            }
        }
        return true;
    }
}

template <typename T>
Array2D<T> Array2D<T>::pow(const T num) {
    /*
    raises the matrix to the given exponent
    
    Input(s):
    num : exponent the matrix is being raised to

    Output(s):
    exp : matrix raised to the given exponent
    */
    Array2D<T> exp = Array2D<T>(arr); // get copy of given array
    for (size_t i = 0; i < this->rows; i++) {
        for (size_t j = 0; j < this->cols; j++) {
            exp.set(pow(exp.get(i,j), num), i, j);
        }
    }
    return exp;
}

template <typename T>
void Array2D<T>::copy_pow(const T num) {
    /*
    raises the matrix to the given exponent
    
    Input(s):
    num : exponent the matrix is being raised to

    Output(s):
    exp : matrix raised to the given exponent
    */
    for (size_t i = 0; i < this->rows; i++) {
        for (size_t j = 0; j < this->cols; j++) {
            this->set(pow(this->get(i,j), num), i, j);
        }
    }
}

template <typename T>
bool Array2D<T>::swap_dim() {
    /*
    swaps a 1D vector from a row to column vector or vice-versa. does nothing
    for non 1D matrices
    
    Input(s): None

    Output(s):
    vec : whether or not the array is a 1D vector
    */
    if ((this->rows == 0) || (this->cols == 0)) {return false;} // deal with non 1-d cases
    else if ((this->rows != 1) && (this->cols != 1)) {return false;}
    size_t rows_temp = this->rows;
    this->rows = this->cols;
    this->cols = rows_temp;
    return true;
}

template <typename T>
void Array2D<T>::print() {
    // prints out array
    for (size_t i = 0; i < this->rows; i++) {
        for (size_t j = 0; j < this->cols; j++) {
            cout << this->get(i,j) << " ";
        }
        cout << endl;
    }
}

template <typename T>
Array2D<T>::~Array2D() {
    // destructor for 2D array
    delete [] this->arr;
}