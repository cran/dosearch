#include <cstdlib>
#include <vector>

#ifndef SET_H
#define SET_H

//#define MAX_SIZE 30
const int MAX_SIZE = 30;

int set_size( const int &set ); /* {
    int size = 0;
    for ( int z=1; z <= MAX_SIZE; z++) {
        if (  set & ( 1 << (z - 1) ) ) size++;
    }
    return size;
}*/

bool in_set( const int &x, const int &set); /* {
    return ( set & ( 1 << (x - 1) ) );
}*/

int full_set(const int &n );

int unary(const int & x);


int set_union(const int & set1, const int & set2);

bool is_subset( const int & set, const int & subset);

// get the next smallest element from a set following previous element
int get_next_element(const int &set, const int &previous);

// get non-empty subsets of a set of size n in order of cardinality (banker's sequence)
std::vector<int> get_subsets(const int &n);

// recursive function to generate banker's sequence
void generate(std::vector<int>& sets, const int& n, int z, int j, int a, const int& b);

#endif
