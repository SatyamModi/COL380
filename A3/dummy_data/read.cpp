#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

void get_binary_0(const char* in_filename, const char* out_filename)
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

void get_binary_1(const char* in_filename, const char* out_filename)
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
 	
	numfs.open("D1.dat", ios::binary | ios::out);
 	
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

void get_binary_2(const char* in_filename, const char* out_filename)
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
 	numfs.open("D2.dat", ios::binary | ios::out);
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

void data_process()
{
	get_binary_0("index.txt", "bin_index.dat");		//done
	get_binary_0("indptr.txt", "bin_ptr.dat");		//done	
	get_binary_0("level_offset.txt", "bin_offset.dat");	
	get_binary_0("level.txt", "bin_level.dat");
	get_binary_1("user.txt", "bin_user.dat");
	get_binary_2("vect.txt", "bin_vect.dat");
}

int main()
{
	// get_binary_user("user.txt", "user_bin.dat");
	data_process();
	return 0;
}
