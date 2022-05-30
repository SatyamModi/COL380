#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>
#include <sys/stat.h>

using namespace std;


bool IsPathExist(string &s)
{
  struct stat buffer;
  return (stat (s.c_str(), &buffer) == 0);
}


void get_binary_0(const char* in_filename,const char* out_filename)
{
	FILE *fs; 
	char ch, buffer[32]; 
	int i = 0, num, j = 0; 
	
	// Openning the file with file handler as fs 
	fs = fopen(in_filename, "r"); 
 	fstream out(out_filename, ios::binary | ios::out);
	// Read the file unless the file encounters an EOF 
	while(1){ 
		// Reads the character where the seeker is currently 
		ch = fgetc(fs); 
		if(ch == EOF){ 
			num = atoi(buffer);
			out.write((char*)&num, 4);
			j++;
			break; 
		} 
		else if(ch == ' ' || ch == '\n'){ 
			num = atoi(buffer);
			out.write((char*)&num, 4); 
			j++; 
			bzero(buffer, 32); 
			i = 0; 
 
			continue; 
		} 
		else{ 
			buffer[i] = ch; 
			i++; 
		} 
	} 
	fclose(fs);
	out.close();
}

void get_binary_1(const char* in_filename, const char* out_filename, string out_path)
{
	FILE *fs; 
	fstream numfs;
 	bool done = false;
 	bool ischar = false;

	char ch, buffer[32]; 
	int i = 0, size = 0, D = 0; 
	float num;
	
	// Openning the file with file handler as fs 
	fs = fopen(in_filename, "r"); 
 	fstream out(out_filename, ios::binary | ios::out);
 	
	numfs.open(out_path + "D1.dat", ios::binary | ios::out);
 	
	// Read the file unless the file encounters an EOF 
	while(1)
	{ 
		// Reads the character where the seeker is currently 
		ch = fgetc(fs); 
		if (ch == EOF || ch == ' ' || ch == '\n')
		{
			if(ch == EOF){ 
				if (ischar)
				{
					num = atof(buffer);
					out.write((char*)&num, 4);
					size++;
					D++;
				}
				numfs.write((char*)&D, 4);
			} 
			//ch == ' ' || ch == '\n'
			else
			{ 
				if (ischar)
				{
					num = atof(buffer);
					out.write((char*)&num, 4); 
					size++; 
					if (ch == '\n')	D++;
					bzero(buffer, 32);
					ischar = false;
				} 
				i = 0; 
			}

			if ((ch == EOF || ch == '\n') && !done)
			{	//get the no. of attrs of a user
				numfs.write((char*)&size, 4);
				done = true;
			}
			if (ch == EOF) break;
		} 
		else
		{ 
			ischar = true;
			buffer[i] = ch; 
			i++; 
		} 
	} 
	fclose(fs);
	out.close();
	numfs.close();
}

void get_binary_2(const char* in_filename, const char* out_filename, string out_path)
{
	FILE *fs; 
 	bool done = false;
 	fstream numfs;

	char ch, buffer[32]; 
	int i = 0, size = 0, D = 0; 
	float num;
	
	// Openning the file with file handler as fs 
	fs = fopen(in_filename, "r"); 
 	fstream out(out_filename, ios::binary | ios::out);
 	numfs.open(out_path + "D2.dat", ios::binary | ios::out);
	// Read the file unless the file encounters an EOF 
	while(1)
	{ 
		// Reads the character where the seeker is currently 
		ch = fgetc(fs); 
		if (ch == EOF || ch == ' ' || ch == '\n')
		{
			if(ch == EOF){ 
				if (done)
				{
					num = atof(buffer);
					out.write((char*)&num, 4);
					D++;
				}
				numfs.write((char*)&D, 4);
				break;
			} 
			//ch == ' ' || ch == '\n'
			else
			{   
				if (done)
				{
					num = atof(buffer);
					out.write((char*)&num, 4); 
					if (ch == '\n') D++;
					bzero(buffer, 32); 
					done = false;
				}
				i = 0; 
			}
		} 
		else
		{ 
			done = true;
			buffer[i] = ch; 
			i++; 
		} 
	} 
	fclose(fs);
	out.close();
	numfs.close();
}

void data_process(string in_path, string out_path)
{
	string in = in_path + "index.txt";
	string out = out_path + "bin_index.dat";
	get_binary_0(in.c_str(), out.c_str());		//done

	in = in_path + "indptr.txt";
	out = out_path + "bin_ptr.dat";

	get_binary_0(in.c_str(), out.c_str());		//done

	in = in_path + "level_offset.txt";
	out = out_path + "bin_offset.dat";	
	get_binary_0(in.c_str(), out.c_str());

	in = in_path + "level.txt";
	out = out_path + "bin_level.dat";	
	get_binary_0(in.c_str(), out.c_str());
	
	in = in_path + "max_level.txt";
	out = out_path + "bin_max_level.dat";	
	get_binary_0(in.c_str(), out.c_str());

	in = in_path + "ep.txt";
	out = out_path + "bin_ep.dat";	
	get_binary_0(in.c_str(), out.c_str());

	in = in_path + "user.txt";
	out = out_path + "bin_user.dat";	
	get_binary_1(in.c_str(), out.c_str(), out_path);

	in = in_path + "vect.txt";
	out = out_path + "bin_vect.dat";	
	get_binary_2(in.c_str(), out.c_str(), out_path);
}

int main(int argc, char* argv[])
{
	// get_binary_user("user.txt", "user_bin.dat");
	string in_path = argv[1];
	string out_path = argv[2];
	
	if (!IsPathExist(out_path)) 
	{
		string cmd = "mkdir " + out_path;
		system(cmd.c_str());
	}
	data_process(in_path, out_path);

	return 0;
}

// with open(f1) as fs:
//     pred = fs.readlines()
// with open(f2) as fs:
//     gh = fs.readlines()

// def makeints(pred):
//     k = []
//     for p in pred:
//         l = str.split(p)
//         print(l)
//         ks = []
//         for v in l:
//             ks.append(int(v))
//         k.append(ks)
//     return k

// pred = makeints(pred)
// gh = makeints(gh)

// precision = 0
// recall = 0

// assert(len(pred) == len(gh))

// for i in range(len(pred)):
//     inter =0
//     for p in pred[i]:
//         if p in gh[i]:
//             inter+=1
//             print("common :", p)
//     if(len(pred[i]) > 0):
//         precision += inter/len(pred[i])
//     if(len(gh[i]) > 0):
//         recall += inter/len(gh[i])

// print("precision ", precision * 100 /len(pred))
// print("recall ", recall * 100/len(gh))