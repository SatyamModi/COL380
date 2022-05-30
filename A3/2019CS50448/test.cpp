#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <chrono>

using namespace std;

struct Comparator 
{
	bool operator()(tuple<int,float> const& t1, tuple<int,float> const& t2)
	{
		return (get<1>(t1) > get<1>(t2));
	}
};

typedef priority_queue<tuple<int,float>, vector<tuple<int,float>>, Comparator> minHeap;
typedef vector<int> _1D_vec_i;
typedef vector<float> _1D_vec_f;
typedef vector<_1D_vec_f> _2D_vec_f;

float cosine_dist(_1D_vec_f& v1, _1D_vec_f& v2) {

	float mod_v1 = 0, mod_v2 = 0;
	float dot_prod = 0;
	int n = v1.size();
	
	for(int i=0; i<n; i++) {
	
		mod_v1 += v1[i]*v1[i];
		mod_v2 += v2[i]*v2[i];
		dot_prod += v1[i]*v2[i];
		
	}
	return 1 - dot_prod/(sqrt(mod_v1)*sqrt(mod_v2));
}

tuple<int,float> max_value_queue(minHeap& pq) {
	
	minHeap pq_new =  pq;
	
	while(pq_new.size() > 1) {
		pq_new.pop();
	}
	
	return pq_new.top();
	
}

minHeap trim_queue(minHeap& pq, int k) 
{
	minHeap pq_new;
	
	while(pq_new.size() != k) {
	
		tuple<int,float> val = pq.top();
		pq_new.push(val);
		pq.pop();
	}
	return pq_new;
}

void show_minHeap(minHeap& pq) 
{
	minHeap pq_new = pq;
	
	while(!pq_new.empty()) {
	
		tuple<int,float> val = pq_new.top();
		
		cout << "(" << get<0>(val) << "," << get<1>(val) << ") ";
		
		pq_new.pop();
	}
	cout << "\n";
}


void vect_transfer(int num_vec, int num_attrs, int size, int rank, _2D_vec_f &vect, string path)
{
	MPI_File in;
	MPI_Status status;
	MPI_Offset offset;
	int i;
    
    int recvcounts[size];
    int displs[size];

    int perNode = num_vec/size;
    offset = perNode*num_attrs*rank;
    int count = perNode*num_attrs;
    for (i = 0; i < size; i++) recvcounts[i] = count;
    for (i = 0; i < size; i++) displs[i] = i*count;

	int count1 = num_vec*num_attrs - (size-1)*(perNode)*num_attrs;
	recvcounts[size-1] = count1;
	if (rank == size-1) count = count1;

    float *b;
    b = new float[num_vec*num_attrs];
    float *a;
    a = new float[count];

    MPI_File_open(MPI_COMM_WORLD, (path + "bin_vect.dat").c_str(), MPI_MODE_RDWR|MPI_MODE_CREATE, MPI_INFO_NULL, &in);
    MPI_File_set_view(in, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL ) ; 

    MPI_File_read_at_all(in, offset, a, count, MPI_FLOAT, &status);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&in);
    
    MPI_Allgatherv(a, count, MPI_FLOAT, b, recvcounts, displs, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
	
	_1D_vec_f temp;

	float* begin = b;
	float* end = b + num_vec*num_attrs;
	
	while (begin != end)
	{
		vect.push_back(_1D_vec_f(begin, begin+num_attrs));
		begin += num_attrs;
	}
}

void user_transfer(int rank, int size, int num_users, int num_attrs, int &k, int &ep, int &max_level, _1D_vec_i &index, _1D_vec_i &indptr, _1D_vec_i &level, _1D_vec_i &level_offset, _2D_vec_f &user, string path)
{
	ifstream file;
	
	//get the entry of ep.txt
	file.open(path + "bin_ep.dat");
	file.read((char*)&ep, 4);
	file.close();
	
	//get the entries of index binary file
	file.open(path + "bin_index.dat");
	int temp_i;
	while(file.read((char *)&temp_i, 4)) 
	{
		index.push_back(temp_i);	
	}
	file.close();

	//get the entries of the indptr binary file
	file.open(path + "bin_ptr.dat");
	while(file.read((char *)&temp_i, 4)) 
	{
		indptr.push_back(temp_i);	
	}
	file.close();

	//get the entries of level binary file
	file.open(path + "bin_level.dat");
	while(file.read((char *)&temp_i, 4)) 
	{
		level.push_back(temp_i);	
	}
	file.close();

	//get the entries of level_offset binary file
	file.open(path + "bin_offset.dat");
	while(file.read((char *)&temp_i, 4)) 
	{
		level_offset.push_back(temp_i);	
	}
	file.close();

	//get the entry of max_level
	file.open(path + "bin_max_level.dat");
	file.read((char*)&max_level, 4);
	file.close();

    int start = rank*(num_users/size)*num_attrs;
    int end = (rank+1)*(num_users/size)*num_attrs;

    if (rank == size-1)
    {
		end = num_users*num_attrs;  	
    }
    
    int count = end-start;
    MPI_File in;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, (path + "bin_user.dat").c_str(), MPI_MODE_RDWR|MPI_MODE_CREATE, MPI_INFO_NULL, &in);
    MPI_File_set_view(in, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL ) ; 

    float *a;
    a = new float[count];
    
    MPI_File_read_at(in, start, a, count, MPI_FLOAT, &status);

    float* begin_ = a;
	float* end_ = a + count;
	
	while (begin_ != end_)
	{
		user.push_back(_1D_vec_f(begin_, begin_+num_attrs));
		begin_ += num_attrs;
	}
}

int **alloc_2d_int(int rows, int cols) {
    int *data = (int *)malloc(rows*cols*sizeof(int));
    int **array= (int **)malloc(rows*sizeof(int*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}
minHeap SearchLayer(_1D_vec_f& q, int k, minHeap& candidates, _1D_vec_i& indptr, _1D_vec_i& index, _1D_vec_i& level_offset, int lc, unordered_set<int>& visited, _2D_vec_f& vect) {

	minHeap top_k = candidates;
	tuple <int,float> top;
	int ep;
	int start;
	int end;
	float _dist;
	
	while(!candidates.empty()) {
	
		top = candidates.top();
		ep = get<0>(top);
		candidates.pop();
		
		start = indptr[ep] + level_offset[lc];
		end = indptr[ep] + level_offset[lc+1];

		
		int px;
		
		for(int i=start; i<end; i++) 
		{
			px = index[i];
			if(visited.find(px) != visited.end() || px == -1) continue;

			visited.insert(px); 
			_dist = cosine_dist(q, vect[px]);
			if(_dist > get<1>(max_value_queue(top_k)) && top_k.size() >= k) continue;
			
			top_k.push(make_tuple(px, _dist));
			
		    if (top_k.size() > k)
				top_k = trim_queue(top_k, k);
		   	
			candidates.push(make_tuple(px, _dist));
		}
	}	
	return top_k;
}

void QueryHNSW(int rank, _1D_vec_f& q, int k, int ep, _1D_vec_i& indptr, _1D_vec_i& index, _1D_vec_i& level_offset , int max_level, _2D_vec_f& vect, int *result) 
{
	minHeap top_k;
	top_k.push(make_tuple(ep, cosine_dist(q, vect[ep])));

	unordered_set<int> visited;
	visited.insert(ep);
	for(int lvl = max_level-1; lvl>=0; lvl--) 
	{
		top_k = SearchLayer(q, k, top_k, indptr, index, level_offset, lvl, visited, vect);
	}

	for (int i = 0; i < k; i++)
	{
		result[i] = get<0>(top_k.top());
		top_k.pop();
	}
}


int main(int argc, char* argv[])
{
	auto begin = std::chrono::high_resolution_clock::now();
	string path = string(argv[1]);
	int k = stoi(argv[2]);
	string output_filename = string(argv[4]);
    int rank, size;

	MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) cout << "No. of nodes: " << size << endl;

    int num_vec, num_users, num_attrs;
    ifstream file;

	file.open(path + "D1.dat");
	file.read((char*)&num_attrs, 4);
    file.read((char*)&num_users, 4);
    file.close();

    file.open(path + "D2.dat");
    file.read((char*)&num_vec, 4);
    file.close();

    //transfer the vect file to each node
    _2D_vec_f vect;
    vect_transfer(num_vec, num_attrs, size, rank, vect, path);

    //transfer the user data to each node
    int ep, max_level;
	_1D_vec_i index, indptr, level, level_offset;
	_2D_vec_f user;
	user_transfer(rank, size, num_users, num_attrs, k, ep, max_level, index, indptr, level, level_offset, user, path);

	//parallelizing across threads in a node
	
	int q_count = user.size();
	int** result;
    result = alloc_2d_int(q_count, k);

	#pragma omp parallel 
	{
		int tid = omp_get_thread_num();
		if (tid == 0 && rank == 0) cout << "num_threads: " << omp_get_num_threads() << endl;
		#pragma omp for
		for (int i = 0; i < q_count; i++)
		{
			QueryHNSW(rank, user[i], k, ep, indptr, index, level_offset , max_level, vect, result[i]);
		}
	}

	#pragma omp barrier

	//transfer the result to node 0
	int **A;
	if (rank == 0)	A = alloc_2d_int(num_users,k);
	
	int recvcounts[size];
	int displs[size];
	
	for (int i = 0; i < size; i++) recvcounts[i] = (num_users/size)*k;
	for (int i = 0; i < size; i++) displs[i] = (num_users/size)*k*i;
	recvcounts[size-1] = num_users*k - (size-1)*(num_users/size)*k;

	if (rank == 0) MPI_Gatherv(&(result[0][0]), recvcounts[rank], MPI_INT, &(A[0][0]), recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	else MPI_Gatherv(&(result[0][0]), recvcounts[rank], MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		ofstream myfile(output_filename);
		for (int i = 0; i < num_users; i++)
		{
			for (int j = 0; j < k; j++)
			{
				// cout << A[i][j] << endl;
				myfile << A[i][j] << " ";
			}
			myfile << "\n";
		}
		myfile.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);

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
