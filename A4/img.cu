#include <bits/stdc++.h>
using namespace::std;

#define root_2 1.414213562373095048801688724209

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__device__ int *data_img;

__device__ int *query_img;

__device__ float *candidate_img;

class pixel {

	public:

	int i,j;
	int r,g,b;
	
	__device__  pixel(int row, int col) {
	
		i = row;
		j = col;
		
	}
	
	__device__ void setRGB(int *img, int m, int n) {
	
		r = img[i*n*3 + j*3 + 0];
		g = img[i*n*3 + j*3 + 1];
		b = img[i*n*3 + j*3 + 2];
		
	}
	
};

struct matchVal {

	float rmsd;
	int x,y;
	int angle;
	
	matchVal(float rmsd_val, int x_val, int y_val, int angle_val) {
	
		rmsd = rmsd_val;
		x = x_val;
		y = y_val;
		angle = angle_val;
	}

};

struct CompareVal {
    bool operator()(matchVal const& p1, matchVal const& p2)
    {
        return p1.rmsd < p2.rmsd;
    }
};


__device__  float RMSD(int m, int n, int i0, int j0, int search_m, int search_n, int *dataImg, int *queryImg) {
		
	//here dataImg and queryImg are different from global data_img and query_img respectively
	// m,n - dimensions of dataImg
	// search_m, search_n - dimensions of queryImg == (search window)
	// lc_x, lc_y - left corner row, col in dataImg
	
	float total = 0;

	for(int i=0; i<search_m; i++) 
	{
		for(int j=0; j<search_n; j++) 
		{
			for(int k=0; k<3; k++) 
			{
				float diff = dataImg[(i+i0)*n*3 + (j+j0)*3 + k] - queryImg[i*search_n*3 + j*3 + k];
				total += diff*diff;
			}
		
		}
	
	}
	
	return sqrt(total/(search_m*search_n*3));
}

__device__ float calcGrayValue_dev(int m, int n, int lx, int ly, int box_m, int box_n) 
{
	// no need to pass img as argument as img is always data_img
	// lx, ly - left corner row and col
	// m, n - dimensions of the image
	// box_m, box_n - dimensions of search box

	float total = 0;
	int count = 0;

	for(int i=lx; i< lx+box_m; i++) 
	{
		for(int j=ly; j< ly+box_n; j++) 
		{
	 		for(int k=0; k<3; k++) 
	 		{
	 		 	total += data_img[i*n*3 + j*3 + k];
	 		 	count++;
	 		}	
		}
	}
	return total/count;
}

__host__ float calcGrayValue_host(int *img, int m, int n, int lx, int ly, int box_m, int box_n) 
{

	// lx, ly - left corner row and col
	// m, n - dimensions of the image
	// box_m, box_n - dimensions of search box

    float total = 0;
    int count = 0;

    for(int i=lx; i< lx+box_m; i++) 
    {
    	for(int j=ly; j< ly+box_n; j++) 
    	{
       		for(int k=0; k<3; k++) 
       		{
       		 	total += img[i*n*3 + j*3 + k];
       		 	count++;
       		}	
    	}
    }

    return total/count;
}

__global__ void kernel0(int i_start, int i_end, int j_start, int j_end, int m1, int n1, int m2, int n2, float th1, float th2, float gray_val_of_query)
{
	
   int tid = blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   int i = tid/n1;
   int j = tid - n1*(tid/n1);
   candidate_img[2*i*n1+2*j] = -1;
	candidate_img[2*i*n1+2*j+1] = -1;	

   if (i >= i_start && i < i_end && j >= j_start && j < j_end)
   {
		float gray_val_of_box = calcGrayValue_dev(m1, n1, i, j, m2, n2);

		if (abs(gray_val_of_query - gray_val_of_box) <= th2)
		{
			float rmsd = RMSD(m1, n1, i, j, m2, n2, data_img, query_img);
			if(rmsd <= th1) 
			{
				// printf("Vertical, Found at lc = (%d,%d) with dist = %f\n", i, j, rmsd);
				candidate_img[2*i*n1+2*j] = rmsd;
				candidate_img[2*i*n1+2*j+1] = 0;
			}
		}
	} 
}

__global__ void kernel1(int i_start, int i_end, int j_start, int j_end, int m1, int n1, int m2, int n2, float th1, float th2, float gray_val_of_query)
{

	int tid = blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   int i = tid/n1;
   int j = tid - n1*(tid/n1);
   candidate_img[2*i*n1+2*j] = -1;
	candidate_img[2*i*n1+2*j+1] = -1;	

   //i and j are the indexes of the lowest point
   if (i >= i_start && i < i_end && j >= j_start && j < j_end)
   {
		int lx = i;
		int ly = j-floor(m2/root_2);
		int box_m = floor((m2+n2)/root_2)+1;
		int box_n = floor(n2/root_2) + floor(m2/root_2)+1;

		float gray_val_of_box = calcGrayValue_dev(m1, n1, lx, ly, box_m, box_n);
		
		if (abs(gray_val_of_query - gray_val_of_box) <= th2)
		{
			float rmsd = 0;
			// printf("i: %d, j: %d, th1: %f, th2: %f, diff: %f\n",i, j, th1, th2, abs(gray_val_of_query - gray_val_of_box) );

			for (int i1 = 0; i1 < m2; i1++)
			{
				for (int j1 = 0; j1 < n2; j1++)
				{
					float x = i + i1/root_2 + j1/root_2;
					float y = j - i1/root_2 + j1/root_2;

					pixel A(floor(x), floor(y));
					pixel B(A.i, A.j+1);
					pixel C(A.i+1, A.j+1);
					pixel D(A.i+1, A.j);
					
					A.setRGB(data_img,m1,n1);
					B.setRGB(data_img,m1,n1);
					C.setRGB(data_img,m1,n1);
					D.setRGB(data_img,m1,n1);
					
					pixel res(-1,-1);
					
					float delta_y = x - floor(x);
					float delta_x = y - floor(y);
					
					res.r = A.r*(1-delta_x)*(1-delta_y) + B.r*delta_x*(1-delta_y) + D.r*(1-delta_x)*delta_y + C.r*delta_x*delta_y;
					res.g = A.g*(1-delta_x)*(1-delta_y) + B.g*delta_x*(1-delta_y) + D.g*(1-delta_x)*delta_y + C.g*delta_x*delta_y;
					res.b = A.b*(1-delta_x)*(1-delta_y) + B.b*delta_x*(1-delta_y) + D.b*(1-delta_x)*delta_y + C.b*delta_x*delta_y;
					
					float delta_r = (query_img[i1*n2*3+j1*3+0]-res.r);
					float delta_g = (query_img[i1*n2*3+j1*3+1]-res.g);
					float delta_b = (query_img[i1*n2*3+j1*3+2]-res.b);

					rmsd += delta_r*delta_r + delta_g*delta_g + delta_b*delta_b;
				}
			}

			rmsd = sqrt(rmsd/(m2*n2*3));
			if (rmsd <= th1) 
			{
				// printf("Angle 45, Found at lc = (%d,%d) with dist = %f\n", i, j, rmsd);	
				candidate_img[2*i*n1+2*j] = rmsd;
				candidate_img[2*i*n1+2*j+1] = 45;
			}
		}
	}
}

__global__ void kernel2(int i_start, int i_end, int j_start, int j_end, int m1, int n1, int m2, int n2, float th1, float th2, float gray_val_of_query)
{
	int tid = blockIdx.x*blockDim.x*blockDim.y + threadIdx.y*blockDim.x + threadIdx.x;
   int i = tid/n1;
   int j = tid - n1*(tid/n1);
   candidate_img[2*i*n1+2*j] = -1;
	candidate_img[2*i*n1+2*j+1] = -1;	

   //i and j are the indexes of the lowest point
   if (i >= i_start && i < i_end && j >= j_start && j < j_end)
   {
		int lx = i - floor(n2/root_2);
		int ly = j;
		int box_m = floor(n2/root_2) + floor(m2/root_2)+1;
		int box_n = floor((m2+n2)/root_2)+1;
		float gray_val_of_box = calcGrayValue_dev(m1, n1, lx, ly, box_m, box_n);
		
		if (abs(gray_val_of_query - gray_val_of_box) <= th2)
		{
			float rmsd = 0;
			// printf("i: %d, j: %d, th1: %f, th2: %f, diff: %f\n",i, j, th1, th2, abs(gray_val_of_query - gray_val_of_box) );

			for (int i1 = 0; i1 < m2; i1++)
			{
				for (int j1 = 0; j1 < n2; j1++)
				{
					float x = i + i1/root_2 - j1/root_2;
					float y = j + i1/root_2 + j1/root_2;

					pixel A(floor(x), floor(y));
					pixel B(A.i, A.j+1);
					pixel C(A.i+1, A.j+1);
					pixel D(A.i+1, A.j);
					
					A.setRGB(data_img,m1,n1);
					B.setRGB(data_img,m1,n1);
					C.setRGB(data_img,m1,n1);
					D.setRGB(data_img,m1,n1);
					
					pixel res(-1,-1);
					
					float delta_y = x - floor(x);
					float delta_x = y - floor(y);
					
					res.r = A.r*(1-delta_x)*(1-delta_y) + B.r*delta_x*(1-delta_y) + D.r*(1-delta_x)*delta_y + C.r*delta_x*delta_y;
					res.g = A.g*(1-delta_x)*(1-delta_y) + B.g*delta_x*(1-delta_y) + D.g*(1-delta_x)*delta_y + C.g*delta_x*delta_y;
					res.b = A.b*(1-delta_x)*(1-delta_y) + B.b*delta_x*(1-delta_y) + D.b*(1-delta_x)*delta_y + C.b*delta_x*delta_y;
					
					float delta_r = (query_img[i1*n2*3+j1*3+0]-res.r);
					float delta_g = (query_img[i1*n2*3+j1*3+1]-res.g);
					float delta_b = (query_img[i1*n2*3+j1*3+2]-res.b);

					rmsd += delta_r*delta_r + delta_g*delta_g + delta_b*delta_b;
				}
			}

			rmsd = sqrt(rmsd/(m2*n2*3));
			if (rmsd <= th1) 
			{
				// printf("Angle -45, Found at lc = (%d,%d) with dist = %f\n", i, j, rmsd);
				candidate_img[2*i*n1+2*j] = rmsd;
				candidate_img[2*i*n1+2*j+1] = -45;		
			}
		}
	}
}

//case when image is +45 rotated
__host__ void searchRotated1(int *dataImg, int m1, int n1, int *queryImg, int m2, int n2, float th1, float th2)
{

	float gray_val_of_query = calcGrayValue_host(queryImg, m2, n2, 0, 0, m2, n2);

	int i_start = 0;
	int i_end = m1 - floor((m2+n2)/root_2) - 1;

	int j_start = ceil(m2/root_2);
	int j_end = n1 - floor(n2/root_2) - 1;

	dim3 dimBlock(32, 32);
	dim3 dimGrid((m1*n1)/1024 + 1, 1);

	kernel1<<<dimGrid , dimBlock>>>(i_start, i_end, j_start, j_end, m1, n1, m2, n2, th1, th2, gray_val_of_query);
	cudaDeviceSynchronize();
}

//case when image is -45 rotated
__host__ void searchRotated2(int *dataImg, int m1, int n1, int *queryImg, int m2, int n2, float th1, float th2)
{
	float gray_val_of_query = calcGrayValue_host(queryImg, m2, n2, 0, 0, m2, n2);

	int i_start = ceil(n2/root_2);
	int i_end = m1 - floor(m2/root_2) - 1;

	int j_start = 0;
	int j_end = n1 - floor((m2+n2)/root_2) - 1;

	dim3 dimBlock(32, 32);
	dim3 dimGrid((m1*n1)/1024 + 1, 1);

	kernel2<<<dimGrid , dimBlock>>>(i_start, i_end, j_start, j_end, m1, n1, m2, n2, th1, th2, gray_val_of_query);
	cudaDeviceSynchronize();
}


__host__ void searchVertical(int *dataImg, int m1, int n1, int *queryImg, int m2, int n2, float th1, float th2) {
	
	//here i and j are lower left corners
	float gray_val_of_query = calcGrayValue_host(queryImg, m2, n2, 0, 0, m2, n2);

	int i_start = 0, i_end = m1-m2+1;
	int j_start = 0, j_end = n1-n2+1;

	dim3 dimBlock(32, 32);
	dim3 dimGrid((m1*n1)/1024 + 1, 1);

	kernel0<<<dimGrid , dimBlock>>>(i_start, i_end, j_start, j_end, m1, n1, m2, n2, th1, th2, gray_val_of_query);
	cudaDeviceSynchronize();

}

__host__ int* readImage(string &filename, int &m, int &n) {

	ifstream file;
	file.open(filename.c_str());

	int x;
	file >> m;
	file >> n;

	int *img = (int*)malloc(3*m*n*sizeof(int));	
	for (int i = 0; i < m; i++)
    {
    	for (int j = 0; j < n; j++)
    	{
    		for (int k = 0; k < 3; k++)
    		{
    			file >> x;
    			img[(m-i-1)*n*3 + j*3 + k] = x;
    		}
    	}
    }

	file.close();
	return img;
}

void getTopMatches(priority_queue<matchVal, vector<matchVal>, CompareVal> &pq, float *candidates, int m1, int n1, int n, bool finalcheck) {

	for(int i=0; i<2*m1*n1; i+=2) 
	{
		
		matchVal currVal(candidates[i], i/(2*n1), (i%(2*n1))/2, (int)candidates[i+1]);

		if(currVal.rmsd < 0) continue;
			
		if(pq.size() < n) 
		{
			pq.push(currVal);	
		}
		else 
		{
			//here topVal has maximum rmsd
			matchVal topVal = pq.top();
				
			if(topVal.rmsd < currVal.rmsd) 
			{
				continue;
			}
			pq.pop();
			pq.push(currVal);
		}
	}

	vector<matchVal> result;
	
	while(!pq.empty() && finalcheck) {
		
		matchVal curr = pq.top();
		//cout << curr.rmsd << " " << curr.x << " " << curr.y << " " << curr.angle << "\n";
		result.push_back(curr);
		pq.pop();
	}
	if (finalcheck)
	{
		ofstream myfile;
		myfile.open("output.txt");
		for (int i = result.size()-1; i >= 0; i--)
		{
			myfile << result[i].x;
			myfile << " ";
			myfile << result[i].y;
			myfile << " ";
			myfile << result[i].angle;
			myfile << "\n";
		}
		myfile.close();
	}
}

int main(int argc, char** argv) 
{
	// th1 => threshold for RMSD
	// th2 => threshold for GrayValue 

	string imagefile = string(argv[1]);
	string queryfile = string(argv[2]);
	float th1=stof(argv[3]), th2=stof(argv[4]);
	int n = stoi(argv[5]);

	int m1,n1;
	int *dataImg = readImage(imagefile, m1, n1);

	int m2, n2;
	int *queryImg = readImage(queryfile, m2, n2);


	int *d_data;
	gpuErrchk(cudaMalloc(&d_data, 3*m1*n1*sizeof(int)));
	gpuErrchk(cudaMemcpy(d_data, &dataImg[0], 3*m1*n1*sizeof(int),cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(data_img, &d_data, sizeof(int*)));

   
	int *d_query;
	gpuErrchk(cudaMalloc(&d_query, 3*m2*n2*sizeof(int)));
	gpuErrchk(cudaMemcpy(d_query, &queryImg[0], 3*m2*n2*sizeof(int),cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(query_img, &d_query, sizeof(int*)));

	float *candidates;
	candidates = (float*)malloc(2*m1*n1*sizeof(float));
	for (int i = 0; i < 2*m1*n1; i++) candidates[i] = -1;

	float *d_candidates;
	gpuErrchk(cudaMalloc(&d_candidates, 2*m1*n1*sizeof(float)));
	gpuErrchk(cudaMemcpy(d_candidates, &candidates[0], 2*m1*n1*sizeof(float),cudaMemcpyHostToDevice));
	gpuErrchk(cudaMemcpyToSymbol(candidate_img, &d_candidates, sizeof(float*)));

	priority_queue<matchVal, vector<matchVal>, CompareVal> pq;

	searchVertical(dataImg, m1, n1, queryImg, m2, n2, th1, th2);
	gpuErrchk(cudaMemcpy(candidates, &d_candidates[0], 2*m1*n1*sizeof(float), cudaMemcpyDeviceToHost));
	getTopMatches(pq, candidates, m1, n1, n, 0);

	searchRotated1(dataImg, m1, n1, queryImg, m2, n2, th1, th2);
	gpuErrchk(cudaMemcpy(candidates, &d_candidates[0], 2*m1*n1*sizeof(float), cudaMemcpyDeviceToHost));
	getTopMatches(pq, candidates, m1, n1, n, 0);

	searchRotated2(dataImg, m1, n1, queryImg, m2, n2, th1, th2);
	gpuErrchk(cudaMemcpy(candidates, &d_candidates[0], 2*m1*n1*sizeof(float), cudaMemcpyDeviceToHost));
	getTopMatches(pq, candidates, m1, n1, n, 1);   

	
	cudaFree(d_query);
	cudaFree(d_candidates);
	cudaFree(d_data);

	return 0;
}
