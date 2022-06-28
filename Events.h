// file for Event struct and related functions

struct Event 
{
    double *times; // pointer to the list of times to trigger the event
    unsigned int num_times; // length of the list
    double period; // how frequently the event should trigger
    double alt; // altitude of the event
};
