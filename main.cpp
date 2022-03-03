//
//  main.cpp
//  matrix-mult
//
//  Created by dwh on 22/10/2021.
//

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <time.h>

const int inf = 3e+8;

typedef struct memo_table {
    bool is_set;
    int num_ops;
    int cut;
} m_table;

typedef struct parenthesis_element {
    int num_l_par_pos;
    int num_r_par_pre;
} p_elem;

m_table ** new_memo_table(int n, int m) {
	m_table ** r = new m_table * [n];
    
    for(int i = 0; i < n; ++i) {
        r[i] = new m_table[m];
    }
    
    return r;
}

void free_memo_table(m_table ** r, int n) {

	for(int i = 0; i < n; ++i)
        delete [] r[i];
    
    delete [] r;
}

int ** int2D(int nx, int ny) {
    int ** p = new int * [nx];

    for(int i = 0; i < nx; ++i)
        p[i] = new int[ny];

    return p;
}

void free_int2D(int ** p, int nx) {

    for(int i = 0; i < nx; ++i)
        delete [] p[i];

    delete [] p;
}

void make_par_tree(int i, int j, m_table ** table, p_elem * par_tree) {
    
    if(j - i > 2) {
        int min_k = table[i][j].cut;
    
        make_par_tree(i, min_k, table, par_tree);
        make_par_tree(min_k, j, table, par_tree);

        if(min_k - i > 1) {
            par_tree[i].num_l_par_pos++;
            par_tree[min_k].num_r_par_pre++;
        }
        
        if(j - min_k > 1) {
            par_tree[min_k].num_l_par_pos++;
            par_tree[j].num_r_par_pre++;
        }
    }
}

int check_num_ops(int i, int j, int p[], m_table ** table) {
    int num_of_ops = 0;
    
    //Case pair of matrices
    if(j - i == 2)
        num_of_ops = p[i] * p[i + 1] * p[j];
    
    //Case product of at least three matrices
    if(j - i > 2) {
        int min_k = table[i][j].cut;
        int num_ops1 = check_num_ops(i, min_k, p, table);
        int num_ops2 = check_num_ops(min_k, j, p, table);
        
        num_of_ops = num_ops1 + num_ops2 + p[i] * p[min_k] * p[j];;
    }
    
    return num_of_ops;
}

int verify_num_ops(int n, int p[], m_table ** table) {
    
    //Verify number of operations by computing matrix product
    int num_ops = check_num_ops(0, n - 1, p, table);
    
    return num_ops;
}

void print_tree(int i, int j, m_table ** table) {
    
    //Initialize parenthesis tree
    int n = j - i + 1;
    p_elem * par_tree = new p_elem[n];
    for(int it = 0; it < n; ++it) {
        par_tree[it].num_l_par_pos = 0;
        par_tree[it].num_r_par_pre = 0;
    }
    
    //Make parenthesis tree
    make_par_tree(i, j, table, par_tree);
    
    //Print tree
    std::cout << "printing matrix product solution:" << std::endl;
    std::cout << "(";
    for(int it = 0; it < n; ++it) {
        int num_l_par_pos = par_tree[it].num_l_par_pos;
        int num_r_par_pre = par_tree[it].num_r_par_pre;
        
        for(int it = 0; it < num_r_par_pre; ++it) { std::cout << ")"; }
        for(int it = 0; it < num_l_par_pos; ++it) { std::cout << "("; }
        if(it < n - 1) { std::cout << "A" << it + 1; }
    }
    std::cout << ")";
    std::cout << std::endl;
    
    //Free data
    delete [] par_tree;
}

int min_ops(int p[], int i, int j, m_table ** table) {
    int min_nops = inf;
    
    //Get value from memo table if available
    if(table[i][j].is_set)
        return table[i][j].num_ops;
    
    //Case single or no matrix
    if(j - i < 2)
        return 0;
    
    //Case product of two matrices
    if(j - i == 2)
        return p[i] * p[i + 1] * p[j];
    
    //Case product of more than two matrices
    int min_k = i;
    if(j - i > 2) {
        for(int k = i + 1; k < j; ++k) {
            int num_ops1 = min_ops(p, i, k, table);
            int num_ops2 = min_ops(p, k, j, table);
            int tot_num_ops = num_ops1 + num_ops2;
            tot_num_ops = tot_num_ops + p[i] * p[k] * p[j];

            if(tot_num_ops < min_nops) {
                min_k = k;
                min_nops = tot_num_ops;
            }
        }
    }
    
    //Set memo table
    table[i][j].is_set = true;
    table[i][j].num_ops = min_nops;
    table[i][j].cut = min_k;

    return min_nops;
}

void print_solution(int n, m_table ** table) {
    
    //Print tree
    print_tree(0, n - 1, table);
}

int minimum_num_ops(int p[], int n, m_table ** table) {
    
    //Initialize memo table
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
        	table[i][j].is_set = false;
        	table[i][j].num_ops = inf;
        }
    }
    
    //Compute minimum number of operations
    int min_num_ops = min_ops(p, 0, n - 1, table);
    
    return min_num_ops;
}

int main(int argc, const char * argv[]) {

    //Input size
    int n = 11;

    //Allocate space for results
    m_table ** table = new_memo_table(n, n);

    //Creating random sequence of matrix dimensions
    int * p = new int[n];
    srand((unsigned) time(NULL));
    for(int i = 0; i < n; ++i) {
        p[i] = rand() % n + 2;
    }

    //Compute minimum number of operations
    int min_num_ops = minimum_num_ops(p, n, table);

    //Verify number of operations
    int min_num_ops_ver = verify_num_ops(n, p, table);

    //Print results
    print_solution(n, table);

    std::cout << "min_num_ops: " << min_num_ops << std::endl;
    std::cout << "min_num_ops_ver: " << min_num_ops_ver << std::endl;

    //Free data
    delete [] p;
    free_memo_table(table, n);

    return 0;
}
