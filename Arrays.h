// 1 and N-d, constant size arrays, along with method for loading vectors

#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
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
        ArrayND(const T val, const array<size_t, N> &dim); // creates array with copies of given value
        ArrayND(const string &filepath); // creates array by loading from .npy file
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
        Array1D(const T val, const size_t len); // creates array with copies of given value
        Array1D(const string &filepath); // load from .npy file
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
ArrayND<T,N>::ArrayND(const string &filepath) {
    /*
    loads an N-d array from a saved .npy type file

    Input(s):
    filepath : string containing relative or absolute path to saved data

    Outpus(s): ArrayND instance
    */

    ifstream data_file; data_file.open(filepath, ios::in | ios::binary);
    if (data_file.is_open()) {
        char waste[8]; // opening bytes of the file that are essentially useless
        unsigned short header_len; // for holding information on the length of the header
        data_file.read(waste, 8); data_file.read((char *) &header_len, 2);
        char * header = new char[header_len]; data_file.read(header, header_len); // get the header
        string header_str = string(header); // convert to C++ string for easier handling
        size_t dim_indx = header_str.find_last_of(':') + 3; // find starting location of dimentions
        size_t next_comma; //  will be pointer for next comma character
        char * num_arr; // charater array that needs to be converted to a number
        size_t num_size; // size of the number
        this->dim = array<size_t,N>(); // dimention array
        size_t curr_dim; // current dimention size
        this->tot_size = 1; // calculate total size of the array
        for (size_t i = 0; i < N; i++) {
            next_comma = header_str.find_first_of(',', dim_indx);
            if ((i == N-1) && (N != 1)) {
                num_size = next_comma - dim_indx - 1; // ignore final bracket
            } else {
                num_size = next_comma - dim_indx;
            }
            num_arr = new char[num_size];
            for (size_t j = 0; j < num_size; j++) { // go through and pull the number
                num_arr[j] = header[j + dim_indx];
            }
            istringstream temp(num_arr); // used to convert to an actual number
            temp >> curr_dim; // convert to size_t
            this->dim[i] = curr_dim; this->tot_size *= curr_dim;
            delete [] num_arr;
        }
        delete [] header;
        
        // drop out if it's a size-zero array
        if (this->tot_size != 0) {
            // otherwise read in the array
            this->arr = new T[this->tot_size]; // create underlying array
            data_file.read((char *) this->arr, (this->tot_size)*sizeof(T)); // read in remaining data
        }
    } else {
        throw invalid_argument(".npy file could not be opened");
    }
}

template <typename T, size_t N>
ArrayND<T,N>::ArrayND(const T val, const array<size_t, N> &dim) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    dim : dimentions of the array

    Output(s): ArrayND instance
    */
    if (N == 0) { // in this case, return zero array
        this->tot_size = 0;
    } else {
        this->tot_size = 1; // calculate total number of elements in the array
        for (size_t i = 0; i < N; i++) {
            this->tot_size *= dim[i];
        }
        this->arr = new T[tot_size]; this->dim = array<size_t,N>(dim);
        // initialize to value
        for (size_t i = 0; i < tot_size; i++) {
            this->arr[i] = val;
        }
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
    size_t skip_mul = this->tot_size;
    for (size_t i = 0; i < N; i++) {
        skip_mul /= (this->dim)[i];
        real_loc += loc[i]*skip_mul;
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
        size_t curr_skip_mul = (this->tot_size)/((this->dim)[0]);
        cout << "[";
        for (size_t i = 0; i < (this->dim)[0]; i++) {
            partial_loc[0] = i;
            this->print_dim(partial_loc, 1, i*curr_skip_mul, curr_skip_mul);
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
        curr_skip_mul /= (this->dim)[num_set];
        for (size_t i = 0; i < (this->dim)[num_set]; i++) {
            partial_loc[num_set] = i;
            this->print_dim(partial_loc, num_set + 1, curr_true_loc + i*curr_skip_mul, curr_skip_mul);
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
Array1D<T>::Array1D(const T val, const size_t len) {
    /*
    initializes an array with all elements set to the given value

    Input(s):
    val : value to fill the array with
    len : length of the array

    Output(s): Array1D instance
    */
    if (len == 0) { // in this case, return zero array
        this->tot_size = 0;
    } else {
        this->arr = new T[len]; this->dim = array<size_t,1>({len}); this->tot_size=len;
        // initialize to zeros
        for (size_t i = 0; i < len; i++) {
            this->arr[i] = val;
        }
    }
}

template <typename T>
Array1D<T>::Array1D(const string &filepath) {
    /*
    loads an 1-d array from a saved .npy type file

    Input(s):
    filepath : string containing relative or absolute path to saved data

    Outpus(s): Array1D instance
    */

    ifstream data_file; data_file.open(filepath, ios::in | ios::binary);
    if (data_file.is_open()) {
        char waste[8]; // opening bytes of the file that are essentially useless
        unsigned short header_len; // for holding information on the length of the header
        data_file.read(waste, 8); data_file.read((char *) &header_len, 2);
        char * header = new char[header_len]; data_file.read(header, header_len); // get the header
        string header_str = string(header); // convert to C++ string for easier handling
        size_t dim_indx = header_str.find_last_of(':') + 3; // find starting location of dimentions
        this->dim = array<size_t,1>(); // dimention array
        size_t curr_dim; // current dimention size
        this->tot_size = 1; // calculate total size of the array
        size_t next_comma = header_str.find_first_of(',', dim_indx); // pointer for next comma character
        size_t num_size = next_comma - dim_indx; // size of the number
        char * num_arr = new char[num_size]; // charater array that needs to be converted to a number
        for (size_t j = 0; j < num_size; j++) { // go through and pull the number
            num_arr[j] = header[j + dim_indx];
        }
        istringstream temp(num_arr); // used to convert to an actual number
        temp >> curr_dim; // convert to size_t
        this->dim[0] = curr_dim; this->tot_size *= curr_dim;
        delete [] num_arr; delete [] header;
        // drop out if it's a size-zero array
        if (this->tot_size != 0) {
            // otherwise read in the array
            this->arr = new T[this->tot_size]; // create underlying array
            data_file.read((char *) this->arr, (this->tot_size)*sizeof(T)); // read in remaining data
        }
    } else {
        throw invalid_argument(".npy file could not be opened");
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

template <typename T>
vector<T> * load_vec(const string &filepath) {
    /*
    loads a vector from a saved .npy type file

    Input(s):
    filepath : string containing relative or absolute path to saved data

    Outpus(s): pointer to vector instance
    */

    ifstream data_file; data_file.open(filepath, ios::in | ios::binary);
    vector<T> * to_return;
    if (data_file.is_open()) {
        char waste[8]; // opening bytes of the file that are essentially useless
        unsigned short header_len; // for holding information on the length of the header
        data_file.read(waste, 8); data_file.read((char *) &header_len, 2);
        char * header = new char[header_len]; data_file.read(header, header_len); // get the header
        string header_str = string(header); // convert to C++ string for easier handling
        size_t dim_indx = header_str.find_last_of(':') + 3; // find starting location of dimentions
        size_t curr_dim; // current dimention size
        size_t tot_size = 1; // calculate total size of the array
        size_t next_comma = header_str.find_first_of(',', dim_indx); // pointer for next comma character
        size_t num_size = next_comma - dim_indx; // size of the number
        char * num_arr = new char[num_size]; // charater array that needs to be converted to a number
        for (size_t j = 0; j < num_size; j++) { // go through and pull the number
            num_arr[j] = header[j + dim_indx];
        }
        istringstream temp(num_arr); // used to convert to an actual number
        temp >> curr_dim; // convert to size_t
        tot_size *= curr_dim;
        delete [] num_arr; delete [] header;
        // drop out if it's a size-zero array
        if (tot_size != 0) {
            // otherwise read in the array
            to_return = new vector<T>(tot_size); // allocate vector
            T * arr = to_return->data(); // access underlying array
            data_file.read((char *) arr, tot_size*sizeof(T)); // read in remaining data
        } else {
            to_return = new vector<T>();
        }
    } else {
        throw invalid_argument(".npy file could not be opened");
    }
    return to_return;
}