#include <string>
#include <mpi.h>
#include <assert.h>
#include "randomizer.hpp"
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <chrono>
using namespace std;

bool sortbysec(const pair<int,int> &a, const pair<int,int> &b)
{
    if (a.second != b.second)
    {
        return a.second > b.second;
    }
    else
    {
        return a.first < b.first;
    }
}

void calculateRank(int rank, int nodePerProcess, vector<vector<int>> adj, int num_nodes, int num_walks, int num_steps, int num_rec, Randomizer r)
{
    fstream file("output.dat", ios::in | ios::out | ios::binary);
    
    int lower = rank*nodePerProcess;
    int upper = min(num_nodes, (rank+1)*nodePerProcess);

    for (int i = lower; i <= upper && i < num_nodes; i++)
    {
        int u = i;

        int outDegree = adj[u].size();
        
        // cout << "Node " << i << "Size: " << outDegree << endl;
        vector<pair<int,int>> COT;
        for (int k = 0; k < num_nodes; k++)
        {
            COT.push_back(make_pair(k, 0));
        }

        //Iterate over all the neighbours of u
        for (int j = 0; j < adj[u].size(); j++)
        {
            int v = adj[u][j];
            int degree = adj[v].size();
            if (degree == 0) continue;

            for (int walk = 0; walk < num_walks; walk++)
            {
                int currNode = v;

                for (int steps = 0; steps < num_steps; steps++)
                {
                    int degree_c = adj[currNode].size();
                    if (degree_c > 0)
                    {
                        int next_step = r.get_random_value(u);
                        // we probabilistically opt to restart
                        if (next_step < 0)
                        {
                            currNode = v;
                        }
                        //go to the next child
                        else
                        {
                            int child = next_step%degree_c;
                            currNode = adj[currNode][child];
                        }
                        COT[currNode].second++;
                    }
                    // currNode has no child, then restart 
                    else
                    {
                        currNode = v;
                    }
                }
            }
        }

        COT[u].second = 0;
        
        for (int j = 0; j < adj[u].size(); j++)
        {
            COT[adj[u][j]].second = 0;
        }
        vector<pair<int,int>> COT_n;
        for (int j = 0;j < COT.size(); j++)
        {
            if (COT[j].second > 0) COT_n.push_back(COT[j]);
        }
        sort(COT_n.begin(), COT_n.end(), sortbysec);

        file.seekp(u*(num_rec*8 + 4));
        outDegree = __builtin_bswap32(outDegree);
        file.write((char*)&outDegree, sizeof(int));
        int count = 0;
        int sz = COT_n.size();
        
        for (int j = 0; j < sz && count < num_rec; j++)
        {
            COT_n[j].first = __builtin_bswap32(COT_n[j].first);
            COT_n[j].second = __builtin_bswap32(COT_n[j].second);
            file.write((char*)&COT_n[j].first, sizeof(int));
            file.write((char*)&COT_n[j].second, sizeof(int));
            count++;
        }
        while (count < num_rec)
        {
            file.write((char*)"NULL", 4);
            file.write((char*)"NULL", 4);
            count++;
        }
    }

    file.close();
}

int main(int argc, char* argv[]){
    assert(argc > 8);
    std::string graph_file = argv[1];
    int num_nodes = std::stoi(argv[2]);
    int num_edges = std::stoi(argv[3]);
    float restart_prob = std::stof(argv[4]);
    int num_steps = std::stoi(argv[5]);
    int num_walks = std::stoi(argv[6]);
    int num_rec = std::stoi(argv[7]);
    int seed = std::stoi(argv[8]);

    Randomizer random_generator(seed, num_nodes, restart_prob);
    int rank, numProcs;
    auto begin = std::chrono::high_resolution_clock::now();
    vector<vector<int>> adj(num_nodes);

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    int i, j;
    int proc;
    int nodePerProcess = num_nodes/numProcs;
    int offset[numProcs];

    //Reading the file
    ifstream file;
    file.open(graph_file);
    while (file.read((char*)&i, 4))
    {
        i = __builtin_bswap32(i);
        file.read((char*)&j, 4);
        j = __builtin_bswap32(j);
        adj[i].push_back(j);
    }
    file.close();
    
    if (rank == 0)
    {
        ofstream out("output.dat", ios::binary | ios::out);
        
        int x = 0;
        for (int i = 0; i < num_nodes; i++)
        {
            x = __builtin_bswap32(x);
            out.write((char*)&x, sizeof(x));
            for (int j = 0; j < num_rec; j++)
            {
                out.write((char*)"NULL", 4);
            }
        }
        out.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // calculateRank(outDegree, adj, procNodes[rank], offset[rank], num_walks, num_steps, num_rec, random_generator);
    calculateRank(rank, nodePerProcess, adj, num_nodes, num_walks, num_steps, num_rec, random_generator);
    MPI_Finalize();
    if (rank == 0)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        double duration = (1e-6 * (std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin)).count());
        cout << "Time taken " <<  duration << "ms" << endl;
    }
    return 0;
}