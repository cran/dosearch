#include "set.h"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <Rcpp.h>

using namespace std;

int set_size( const int &set ) {
    int size = 0;
    for ( int z=1; z <= MAX_SIZE; z++) {
        if (  set & ( 1 << (z - 1) ) ) size++;
    }
    return size;
}

bool in_set( const int &x, const int &set) {
    return ( set & ( 1 << (x - 1) ) );
}

int set_union(const int & set1, const int & set2) {
    return(set1 | set2);
}

int full_set(const int &n ) {
    return((1 << n)-1);
}

int unary(const int & x) {
    return( (1 << (x-1)) );
}

bool is_subset( const int & set, const int & subset) {
    return ((set & subset) == subset);
}

int get_next_element(const int& set, const int& previous) {
    for ( int z = 1 + previous; z <= MAX_SIZE; z++ ) {
        if ( set & ( 1 << (z - 1) ) ) return(z);
    }
    return(0);
}

vector<int> get_subsets(const int &n) {
    vector<int> sets;
    for ( int i = 1; i <= n; i++ ) {
        generate(sets, n, 0, 0, 0, i);
    }
    return sets;
}

void generate(std::vector<int>& sets, const int& n, int z, int j, int a, const int& b) {
    if ( a < b ) {
        if ( a == 0 ) {
            for ( int i = 1; i <= n; i++ ) {
                int zz = 1 << (i - 1);
                generate(sets, n, zz, i, a + 1, b);
            }
        } else {
            for ( int i = j + 1; i <= n; i++ ) {
                int zz = z + (1 << (i - 1));
                generate(sets, n, zz, i, a + 1, b);
            }
        }
    } else {
        sets.push_back(z);
    }
}
