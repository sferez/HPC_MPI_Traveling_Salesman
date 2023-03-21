// B&B PARALLEL (MPI) S392371
/**
 * @file mainMPI.cpp
 * @brief Parallel implementation of the Branch and Bound algorithm using MPI
 * @author Simeon FEREZ S392371
 * @version 1.0
 * @date 2023-02-0
 *
 * */

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <mpi.h>
#include <fstream>
#include <climits>
#include <thread>

std::vector<int>  distances; // 1D vector to store the distances between cities
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
                    for (int y = 0; y < num_cities; y++) {
                        if (y==start || y==i || y==j || y==z || y==x) continue;
                        std::vector<int> path;
                        path.push_back(start);
                        path.push_back(i);
                        path.push_back(j);
                        path.push_back(z);
                        path.push_back(x);
                        path.push_back(y);
                        paths.push_back(path);
                    }
                }
            }
        }
    }
}

/**@brief Function to execute the B&B prunning algorithm to find the shortest path in a deep first search manner
 * @param path: the current path that is being explored
 * @param visited: a vector of booleans to keep track of the cities visited
 * @param curr_distance: the current distance of the explored path so far
 * @param local_min_distance: the local minimum distance found so far
 * @param local_min_path: the local minimum path found so far
 * @return: void
 **/
void dfs(std::vector<int>& path, std::vector<int> &visited, int &curr_distance, int &local_min_distance, std::vector<int> &local_min_path) {

    // Base case: all cities in the assigned range have been visited
    if (path.size() == num_cities) {
        // Add the distance from the last city back to the starting city
        int dist = curr_distance + get_distance(path[path.size() - 1], path.back());
        if (dist < local_min_distance) {
            local_min_distance = dist;
            local_min_path = path;
        }
        return;
    }

    // Try visiting all unvisited cities in the assigned range
    for (int i = 0; i < num_cities; i++) {
        if (!visited[i]) {
            int prev_distance = curr_distance;
            curr_distance += get_distance(path.back(), i);
            // Prune the search tree if the current distance is already greater than the local minimum distance
            if (curr_distance >= local_min_distance) {
                curr_distance = prev_distance;
                continue;
            }
            visited[i] = true;
            path.push_back(i);
            dfs(path,  visited, curr_distance, local_min_distance, local_min_path);
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

    // Initialize MPI
    MPI_Init(NULL, NULL);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Define variables
    int buffer;
    MPI_Request request;
    double start_time = MPI_Wtime();
    int local_min_distance = INT_MAX;
    int global_min_distance;
    std::vector<int> local_min_path;
    std::vector<std::vector<int> > paths;

    num_cities = read_file(filename); // read the file

    // Generate random starting city
    int first_city ;
    if (rank == 0) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, num_cities - 1);
        first_city = dis(gen);
    }
    MPI_Bcast(&first_city, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast the starting city to all processes

    // Create all possible pre-paths to be explored and shuffle them
    create_paths(paths,first_city);
    std::mt19937 g(num_procs);
    std::shuffle(paths.begin(), paths.end(), g);

    // Print information about the program
    if (rank ==0){
        std::cout << "B&B PARALLEL" << std::endl;
        std::cout << "Number of processes: " << num_procs << std::endl;
        std::cout << "Number of pre-paths to explore: " << paths.size() << std::endl;
        std::cout << "Number of cities: " << num_cities << std::endl;
        std::cout << "Starting city: " << first_city+1 << std::endl;
    }

    //divide the pre-paths between the processes
    int chunk = paths.size()/num_procs;
    double idle_time_total= 0;
    int start = rank*chunk;
    int end = (rank+1)*chunk-1;
    int rest = paths.size()%num_procs;
    int iter = 0;
    if(rank == num_procs-1){
        end = paths.size()-1;
    }

    // Explore the pre-paths assigned to the process
    for (int i = start; i <= end; i++) {

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

        // Explore the pre-path
        dfs( path, visited, curr_distance, local_min_distance, local_min_path);

//        // Reduce the local_min_distance to every process using a non-blocking communication
//        if (i < paths.size() - rest-1 && iter%720==72) {
//            buffer = local_min_distance; // store the local_min_distance in a buffer to be sent
//            MPI_Iallreduce(&buffer, &global_min_distance, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &request);
//        }
//
//        // Reduce the local_min_distance to every process using a non-blocking communication
//        if (i < paths.size() - rest-1 && iter%720==360) {
//            // Wait for the non-blocking communication to finish
//            double idle_time = MPI_Wtime();
//            MPI_Wait(&request, MPI_STATUS_IGNORE);
//            idle_time_total += MPI_Wtime() - idle_time;
//            // Update the local_min_distance
//            if (global_min_distance < local_min_distance) {
//                local_min_distance = global_min_distance;
//            }
//        }

        if (i < paths.size() - rest-1 && iter%720==72){
            // Do a blocking communication to update the local_min_distance
            buffer = local_min_distance;
            double idle_time = MPI_Wtime();
            MPI_Allreduce(&buffer, &global_min_distance, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            idle_time_total += MPI_Wtime() - idle_time;
            // Update the local_min_distance
            if (global_min_distance < local_min_distance) {
                local_min_distance = global_min_distance;
            }
        }

        iter++;
    }

    // Recalculate the local_min_distance for the local_min_path
    local_min_distance = 0;
    for (int j = 0; j < num_cities - 1; j++) {
        local_min_distance += get_distance(local_min_path[j], local_min_path[j + 1]);
    }

    // Gather minimum distances and paths from all ranks
    std::vector<int> global_min_path;
    if (rank == 0) {
        global_min_path.resize(num_cities);
    }
    int* min_distances = new int[num_procs];
    MPI_Allgather(&local_min_distance, 1, MPI_INT, min_distances, 1, MPI_INT, MPI_COMM_WORLD);

    // Find rank with global minimum distance
    int global_min_rank = std::distance(min_distances, std::min_element(min_distances, min_distances + num_procs));

    // Receive global minimum path from rank with global minimum distance
    if (rank == 0) {
        if (global_min_rank != 0) {
            MPI_Recv(global_min_path.data(), num_cities, MPI_INT, global_min_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            global_min_path = local_min_path;
        }
        global_min_distance = min_distances[global_min_rank];
    } else if (rank == global_min_rank) {
        MPI_Send(local_min_path.data(), num_cities, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    // Free min_distances array
    delete[] min_distances;

    // Print final results
    if(rank==0){
        double end_time = MPI_Wtime();
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "Time: " << end_time - start_time << std::endl;
        std::cout << "Time x Procs: " << (end_time - start_time)*num_procs << std::endl;
        std::cout << "Minimum distance: " << global_min_distance << std::endl;
        std::cout << "Minimum path: ";
        for (int i = 0; i < global_min_path.size(); i++) {
            std::cout << global_min_path[i]+1 << ", ";
        }
        std::cout << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
    }

    // Print idle time
    MPI_Barrier(MPI_COMM_WORLD);
    std::this_thread::sleep_for(std::chrono::milliseconds(100*rank));
    std::cout << "Rank : " << rank << " ; Idle time: " << idle_time_total << std::endl;
//    std::cout << idle_time_total << std::endl;
    MPI_Finalize();

    return 0;
}
