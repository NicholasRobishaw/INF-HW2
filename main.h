#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <ctime>
#include <memory>

using namespace std;

void initial_Construct();
void query_Constructor(const int new_Size, const bool is_First, const string new_Query);
bool read_Qurey(const string& file_Name);
int searchQuery(bool linear);
void sortFragments(char** query_Data, int size, bool ms);
bool file_reader(const string& file_Name);
void qurey_Deconstructor(char** arr_To_Destroy, int size);
void genome_Deconstructor(char* genome_arr, int size);
void copy_String(char* main_Str, long location, long total_Size, const string to_Add);
char* resize_Char_Arr(char* old_Char_Arr, long old_size, long new_size, bool is_First);
string* resize_Str_Arr(string* old_Str_Arr, long old_size, long new_size, bool is_First);
int* resize_Int_Arr(int* old_Int_Arr, long old_size, long new_size, bool is_First);
void mergeSort( char** arr, int low, int high);
void merge(char** arr, int low, int mid, int high);
void bubbleSort(char** arr, int len);
void swap(int *oneVal, int *otherVal);
long linear_Search(string target, char** query, long start, long end, long query_Length);
long binary_Search(string target, char** query, long low, long high);
void copy_Col(char** query, int row, char fragment[]);
