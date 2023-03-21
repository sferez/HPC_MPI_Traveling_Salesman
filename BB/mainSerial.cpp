// B&B SERIAL S392371

/**
 * @file mainMPI.cpp
 * @brief Serial implementation of the Branch and Bound algorithm using MPI
 * @author Simeon FEREZ S392371
 * @version 1.0
 * @date 2023-02-0
 *
 * */

#include <iostream>
#include <algorithm>
#include <chrono>
#include <vector>
#include <climits>
#include <fstream>
#include <random>

int min_distance = INT_MAX;
std::vector<int> min_path;

std::vector<int> distances; // 1D vector to store the distances between cities
int num_cities; // number of cities

/** @brief: Function to get the distance between two cities
 * @param city1: the first city
 * @param city2: the second city
 * @return: the distance between the two cities
 **/
int get_distance(int city1, int city2) {
    if (city1 >= num_cities || city2 >= num_cities) {
        std::cerr << "Error: Invalid city number." << std::endl;
        return 999;
    }
    return distances[city1 * num_cities + city2];
}

/** @brief Function to read the file and store the distances between cities
 * @param filename: the name of the file
 * @return: the number of cities
 **/
int read_file(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file." << std::endl;
        return 0;
    }
    infile >> num_cities;

    distances.resize(num_cities * num_cities);
    for (int i = 0; i < num_cities; i++) {
        for (int j = 0; j < i; j++) {
            int distance;
            infile >> distance;
            distances[i * num_cities + j] = distance;
            distances[j * num_cities + i] = distance;
        }
    }

    infile.close();
    return num_cities;
}

/** @brief Function to create all possible pre-paths to be explore by the brute force algorithm
 * @param paths: a vector of vectors to store all possible paths
 * @param start: the starting city
 * @return: void
 **/
void create_paths( std::vector<std::vector<int> > &paths, int start) {
    for (int i = 0; i < num_cities; i++) {
        if (i==start) continue;
        for (int j = 0; j < num_cities; j++) {
            if (j==start || j==i) continue;
            for (int z = 0; z < num_cities; z++) {
                if (z==start || z==i || z==j) continue;
                for (int x = 0; x < num_cities; x++) {
                    if (x==start || x==i || x==j || x==z) continue;
                    std::vector<int> path;
                    path.push_back(start);
                    path.push_back(i);
                    path.push_back(j);
                    path.push_back(z);
                    path.push_back(x);
                    paths.push_back(path);
                }
            }
        }
    }
}

/**@brief Function to execute the B&B prunning algorithm to find the shortest path in a deep first search manner
 * @param path: the current path that is being explored
 * @param visited: a vector of booleans to keep track of the cities visited
 * @param curr_distance: the current distance of the explored path so far
 * @return: void
 * */
void dfs(std::vector<int> &path, std::vector<int> &visited, int &curr_distance) {

    // Base case: all cities have been visited
    if (path.size() == num_cities) {
        // Add the distance from the last city back to the starting city
        int dist = curr_distance + get_distance(path[path.size() - 1], path.back());
        if (dist < min_distance) {
            min_distance = dist;
            min_path = path;
        }
        return;
    }

    // Try visiting all unvisited cities
    for (int i = 0; i < num_cities; i++) {
        if (!visited[i]) {
            int prev_distance = curr_distance;
            curr_distance += get_distance(path.back(), i);
            // Prunning: if the current distance is already greater than the minimum distance, then stop exploring this path
            if (curr_distance >= min_distance) {
                curr_distance = prev_distance;
                continue;
            }
            visited[i] = true;
            path.push_back(i);
            dfs(path, visited, curr_distance);
            path.pop_back();
            visited[i] = false;
            curr_distance = prev_distance;
        }
    }
}

int main(int argc, char* argv[]) {

    std::string filename; // name of the file
    if (argc == 2) {
        filename = argv[1];
    }
    else {
        std::cout << "Please enter a filename" << std::endl;
        return 0;
    }

    num_cities = read_file(filename);
    clock_t start = clock();

    // Generate random starting city
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, num_cities - 1);
    int first_city = dis(gen);

    //Print information about the program
    std::cout << "B&B SERIAL" << std::endl;
    std::cout << "Number of cities: " << num_cities << std::endl;
    std::cout << "Starting city: " << first_city +1 << std::endl;

    // Create all possible pre-paths to be explored and shuffle them
    std::vector<std::vector<int> > paths;
    create_paths( paths,first_city);
    std::shuffle(paths.begin(), paths.end(), gen);

    // Explore all possible pre-paths
    for (int i = 0; i < paths.size(); i++) {

        std::vector<int> path = paths[i];
        std::vector<int> visited(num_cities, false);
        int curr_distance=0 ;

        // Initialize the curr-distance and visited vector for the current pre-path
        for (int j = 0; j < path.size(); j++) {
            if (j<path.size()-1) {
                curr_distance += get_distance(path[j], path[j + 1]);
            }
            visited[path[j]] = true;
        }

        // Explore the current pre-path
        dfs(path,  visited, curr_distance);
    }

    clock_t end = clock();
    double elapsed_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

    // Print the final results
    std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
    std::cout << "Minimum distance: " << min_distance << std::endl;
    std::cout << "Minimum path: ";
    for (int i = 0; i < min_path.size(); i++) {
        std::cout << min_path[i]+1 << ", ";
    }
    std::cout << std::endl;

    return 0;
}
