//
// Created by taratt on 2/5/23.
//
/**
 * This library contains functions to perform Computations on a matrix's dependency graph
*/

#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/named_function_params.hpp>

#include <typeinfo>
#include <boost/graph/connected_components.hpp>
#include "queue"
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, boost::no_property> Graph;

/**
 * Implements a recursive approach to visit vertices of a graph using DFS
 *
 * @param reached The already reached (visited) vertices
 * @param visited An array that saves 0 for indices that their corresponding vertex has not been visited yet and 1 for those who has
 * @param new_vertex The vertex that the algorithm starts visiting from
 * @param graph The dependency graph
 */
void visit_vertices_recursive(set<int> &reached, int * visited, int new_vertex, Graph graph){
    reached.insert(new_vertex);
    visited[new_vertex] = 1;

    //Iterating over outgoing edges of the current vertex
    for (auto ed : make_iterator_range(boost::out_edges(new_vertex, graph))){
        if (!visited[ed.m_target]){
            visit_vertices_recursive(reached, visited, ed.m_target, graph);
        }
    }
}
/**
 * Implements an iterative approach to visit vertices of a graph using BFS
 *
 * @param reached The already reached (visited) vertices
 * @param visited An array that saves 0 for indices that their corresponding vertex has not been visited yet and 1 for those who has
 * @param new_vertex The vertex that the algorithm starts visiting from
 * @param graph The dependency graph
 */
void visit_vertices_iterative(set<int> &reached, int * visited, int new_vertex, Graph graph){
    queue<int> q;
    reached.insert(new_vertex);
    visited[new_vertex] = 1;
    q.push(new_vertex);
    while (!q.empty())
    {
        new_vertex = q.front();
        q.pop();
    //Iterating over outgoing edges of the current vertex
        for (auto ed : make_iterator_range(boost::out_edges(new_vertex, graph)))
        {
            if (!visited[ed.m_target])
            {
                visited[ed.m_target] = true;
                reached.insert(ed.m_target);
                q.push(ed.m_target);
            }
        }
    }
}

/**
 * Implements an iterative approach to visit vertices of the graph and determine the dependencies between them using BFS
 * A vertex is depended on another when there is a path from the latter to the first in the directed graph
 *
 * @param dependency_items A vector of sets that saves the vertices that each vertex is depended on
 * @param new_vertex The vertex that the algorithm starts visiting from
 * @param graph The dependency graph
 */
void visit_dependencies_iterative(vector<set<int>> &dependency_items, int new_vertex, Graph graph){
    queue<int> q;

    q.push(new_vertex);
    while (!q.empty())
    {
        new_vertex = q.front();
        q.pop();
    //Iterating over outgoing edges of the current vertex
        for (auto ed : make_iterator_range(boost::out_edges(new_vertex, graph)))
        {
                //Saving the source of the edge in the dependencies of the target
                dependency_items[ed.m_target].insert(ed.m_source);
                q.push(ed.m_target);
        }
    }
}

/**
 * Creates the dependency graph of a matrix (which is a directed graph) where each edge from the vertex i to the vertex j,
 * represents a nonzero element of the matrix in the row j and the column i
 *
 * @param matrix The CSC format of the matrix
 * @return graph The dependency graph
 */

Graph create_dependency_graph(CSC matrix){
    Graph graph;
    //Looping over all nonzero elements of the matrix
    for (int j = 0; j<matrix.dim; j++){
        for (int p = matrix.col_ptr[j]; p<matrix.col_ptr[j+1]; p++) {
            if (j!=matrix.row_ind[p])
                //adding the corresponding edge to the graph
                add_edge(j, matrix.row_ind[p],graph);
        }
    }

    return graph;
}

/**
 * Calculated the columns that take part in the triangular solve algorithm
 *
 * @param matrix The CSC format of the matrix
 * @param x The right hand side vector of the system
 * @return a set of the column numbers that take part in the triangular solve algorithm
 */
set <int> get_reach(CSC matrix, Vector_b x){
    Graph graph = create_dependency_graph(matrix);
    int * visited = new int [matrix.dim]();

    set <int> mattered_columns;

    //All rows that contain a nonzero value in x and columns that can be reached from them take part.
    //Calculating reached columns from nonzero rows of x
    for (int i = 0; i < x.nonzeros; ++i) {
        visit_vertices_iterative(mattered_columns, visited, x.nonzero_rows[i],graph);
    }

    return mattered_columns;
}

vector<vector <int>> create_levels(CSC matrix, Vector_b x){
    set <int> reached_columns = get_reach(matrix, x);
    Graph graph = create_dependency_graph(matrix);
    set<int >::iterator reach_it ;
    int * levels = new int [matrix.dim]();
    vector<int> reached_vec (reached_columns.begin(), reached_columns.end());
    for (reach_it = reached_columns.begin() ; reach_it != reached_columns.end() ; reach_it++ ){
        int j = *reach_it;
        for (auto ed : make_iterator_range(boost::out_edges(j, graph))){
            levels[ed.m_target] = max(levels[j]+1, levels[ed.m_target]);
        }
    }

    vector<vector<int>> partitions;
    partitions.emplace_back();
    for(int i=0; i<reached_columns.size(); i++){
        int x = reached_vec[i];
        while(levels[x]> partitions.size()-1){
            partitions.emplace_back();
        }
        partitions[levels[x]].push_back(x);
    }
    return (partitions);
}

void get_dependencies(CSC matrix, Vector_b x,set <int> &reached_columns, int * &dependency_nums){
    Graph graph = create_dependency_graph(matrix);
    vector<set<int>> dependency_items;
    reached_columns = get_reach(matrix,x);
    vector<int> reached_vec (reached_columns.begin(), reached_columns.end());

    for (int i = 0; i < matrix.dim; ++i) {
        dependency_items.emplace_back();
    }
    for (int i=0 ;i<x.nonzeros ; i++){

        visit_dependencies_iterative(dependency_items, x.nonzero_rows[i],graph);
    }
    for (int i = 0; i < reached_columns.size(); ++i) {
        dependency_nums[reached_vec[i]] = dependency_items[reached_vec[i]].size();
    }


}

