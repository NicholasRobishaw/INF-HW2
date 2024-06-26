// includes
#include "main.h"


using namespace std;
// gloabals


// functions

// Queries_AR class
class Queries_AR {
    public:
    // storage for query dataset (2D Array)
    char** query_Data;
    int query_Size;

    // storage for subject dataset
    char* genome_Data;
    int genome_Size;

    int* found_Frags;
    int frag_Size;
    
    const bool debug_Statements = true;
    const int fragment_Size = 32;
    

    // const int fragment_Size = 32;

   

    // default constructor function
    void initial_Construct(){
        query_Data = nullptr;
        genome_Data = nullptr;

        query_Size = 0;
        genome_Size = 0;

        found_Frags = nullptr;
        frag_Size = 0;
    }

    // custom constructor function for resizing the mulitdimesional array
    void query_Constructor(const int new_Size, const bool is_First, const string new_Query){
        int index, row, col;
        char** new_Temp = nullptr;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "entered query_Constructor function\n";
            cout << "New size: " << new_Size << endl;
        }

        // allocate initial single array of pointers
        new_Temp = new char*[new_Size];

        if(debug_Statements){
            cout << "Array of pointers has been dynamically created\n";
            cout << "Now iterating thorugh rows and allocating 2nd Dimension mem\n";
        }

        // loop thorugh initial array and allocate spaces for each index needed
        for(index=0;index<new_Size;index++){
            // give the pointer some memory
            new_Temp[index] = new char[fragment_Size];
        }
        
        if(debug_Statements){
            cout << "rows have been successfully allocated\n";
        }

        // check if this is NOT the first query
        if(!is_First){
            
            if(debug_Statements){
                cout << "now copying over old array values to new array\n";
            }   
            
            // copy data into new array
            for(row=0;row<new_Size;row++){
                for(col=0;col<fragment_Size;col++){
                    new_Temp[row][col] = query_Data[row][col];
                }
            }


            if(debug_Statements){
                cout << "copy completed, calling destructor\n";
            } 
            
            // delete old array
            qurey_Deconstructor(query_Data, query_Size);

            if(debug_Statements){
                cout << "destructor completed\n";
            } 

        }

        if(debug_Statements){
            cout << "adding new items to end of array\n";
        } 

        // add new data into the new spot
        for(index=0;index<fragment_Size;index++){
            // ary[y*sizeX + x];
            new_Temp[new_Size-1][index] = new_Query[index];
        }


        if(debug_Statements){
            cout << "new items added, setting ptr values\n";
        } 

        // set pointer to new array
        query_Data = new_Temp;
        new_Temp = nullptr;

    }

    // read query dataset file function
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string current_Line;
        int query_Num = 0;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        // error handling for mainly bad_allocs or segmentation faults
        try {
            // iterate through file contents line by line
            while (getline(file, current_Line)) {
                // check if we are at a valid fragment
                if(current_Line[0] != '>' && !current_Line.empty()){
                    // increment query count
                    query_Num++;
                    
                    
                    if(debug_Statements){
                          cout << "rd_Qry-Q_NUM: " << query_Num << endl;
                          cout << "rd-Qry-curLn: " << current_Line << endl;
                    } 
                    

                    // resize query array
                    query_Constructor(query_Num, query_Num == 1, current_Line);
                    
                    if(debug_Statements){
                        cout << "rd_Qry: query that returned\n";
                        for(int dindex = 0; dindex < query_Num; dindex++){
                            for(int dcol = 0; dcol<fragment_Size; dcol++){
                                cout << query_Data[dindex][dcol];
                            }
                            
                            cout << endl;
                        }
                    }

                    // increment query size variable in class
                    query_Size++;
                }

                
            }

            // close the file
            file.close();

            // return sucess
            return true;

        // catch any bad alloc returns and throw them up to main funciton
        } catch (const std::bad_alloc& e){
            // send back error message and return failure
            cerr << "Memory allocation failed: " << e.what() << endl;
            return false;
        
        // catch any other errors that occur and throw them up to main funciton
        } catch(const std::exception& e){
            // send back error message and return failure
            cerr << "Exception occurred: " << e.what() << endl;
            return false;
        }

        return 0;
    }

    // search function
    int searchQuery(bool linear){
        long index = 0, high, target_Index = 0, target_Index_Src = 0, match_Index;
        bool frag_Found = false;
        int query_Index = 0;
        string target = "";

        // iterate through 32-character long fragments of the subject dataset, searching for each one in the query dataset
        while(index < genome_Size && !frag_Found){

            high = index+=fragment_Size;

            if( high < genome_Size ){
                target_Index_Src = 0;
                target = "";

                for(target_Index = index; target_Index<high; target_Index++){
                    target[target_Index_Src] = genome_Data[target_Index];
                    target_Index_Src;
                }


                // linear Search?
                if(linear){
                    match_Index = linear_Search(target, query_Data, index, high, query_Size);
                }
                // Binary Search?
                else{
                    match_Index = binary_Search(target, query_Data, index, high);
                }

                resize_Int_Arr(found_Frags, frag_Size, frag_Size+1, frag_Size == 0);
                frag_Size++;
                found_Frags[frag_Size] = match_Index;

                // copy_String(found_Frags, genome_Index, new_genome_size, temp_Genome )
                query_Index++;
            }



            index+=fragment_Size;
        }

        // return index of match or negative value if no hit was found
        return 0;
    }

    // sort fragments function
    void sortFragments(char** query_Data, int size, bool ms){
        // variables

        // sort all character fragments in alpahbetical order (lexicographic) order (logorithmic runtime!)
        // check for merge sort code
        if(ms){
            mergeSort(query_Data, 0, size);
        }
        // otherwise assume bubble sort
        else{
            bubbleSort(query_Data, size);
        }

    }

    // function to read in the genome data
    // function for reading genome file ( LOOK AT THIS WHEN WE GET THE INPUT FILE )
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string current_Scaffold_Name;
        string current_Line;
        string temp_Genome = "";
        long new_genome_size;
        long protein_Count;
        long genome_Index;
        long index;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }

        // error handling for mainly bad_allocs or segmentation faults
        try {
            // iterate through file contents line by line
            while (getline(file, current_Line)) {
                // check if we hit a header
                if(current_Line[0] == '>'){
                    // increment scaffold count by 1
                    genome_Size++;
                    
                    // check and make sure that there is a genome thats ready to be processed
                    if( !temp_Genome.empty()){
                        // grab the protein count of the current scaffold
                        protein_Count = temp_Genome.length();
                        
                        // calculate new size for genome array
                        genome_Index = genome_Size;
                        new_genome_size = genome_Size + protein_Count;
        
                        // toggle for displaying debug statements
                        if( debug_Statements ){
                            cout << "Resizing genome array to fit current scaffold\n";
                        }
        
                        // resize the genome array to fit last scaffold and copy data into new array
                        genome_Data = resize_Char_Arr(genome_Data, genome_Size, new_genome_size, 
                                                                    genome_Size == 0);
        
                        // toggle for displaying debug statements
                        if( debug_Statements ){
                            cout << "Resize completed adding new values to genome\n";
                        }
        
                        // concatenate the temp_Genome string into the genome character array
                        copy_String(genome_Data, genome_Index, new_genome_size, temp_Genome );
        
                        // set the new tracking size for the genome
                        genome_Size = new_genome_size;
        
                        // reset the temp genome string to empty
                        temp_Genome = "";
                    }
                }

                // check if the line is not empty or a header
                else if (!current_Line.empty() && current_Line[0] != '>') {
                    // add to the current scaffold string
                    temp_Genome += current_Line;
                }
            }

            // post append information of last scaffold
            if (!temp_Genome.empty()) {
                // toggle for displaying debug statements
                if( debug_Statements ){
                    cout << "Last scaffold post loop being processed\n";
                }

                // grab the length of the last scaffold
                protein_Count = temp_Genome.length();
                
                // calculate new genome array size
                new_genome_size = genome_Size + protein_Count;
                
                // toggle for displaying debug statements
                if( debug_Statements ) {
                    cout << "Resizing genome arr for last scaffold\n";
                }

                // resize the genome array to fit last scaffold and copy data into new array
                genome_Data = resize_Char_Arr(genome_Data, genome_Size, new_genome_size, genome_Size == 0);

                // concatenate the temp_Genome string into the genome character array
                copy_String(genome_Data, genome_Index, new_genome_size, temp_Genome );

                // set genome size tracker
                genome_Size = new_genome_size;

            }

            // close the file
            file.close();

            // return sucess
            return true;

        // catch any bad alloc returns and throw them up to main funciton
        } catch (const std::bad_alloc& e){
            // send back error message and return failure
            cerr << "Memory allocation failed: " << e.what() << endl;
            return false;
        
        // catch any other errors that occur and throw them up to main funciton
        } catch(const std::exception& e){
            // send back error message and return failure
            cerr << "Exception occurred: " << e.what() << endl;
            return false;
        }
    }// end of file reader function


    // deconstructor function for query
    void qurey_Deconstructor(char** arr_To_Destroy, int size){
        int index;

        for(index=0;index<size;index++){
            delete[] arr_To_Destroy[index];
        }

        delete[] arr_To_Destroy;
        arr_To_Destroy=nullptr;
    }

    void genome_Deconstructor(char* genome_arr, int size){
        delete[] genome_arr;
        genome_arr=nullptr;
    }

    void copy_String(char* main_Str, long location, long total_Size, const string to_Add){
        int index = 0;
        
        for(location; location < total_Size; location++){
            main_Str[location] = to_Add[index];
            index++;
        }
    }

    char* resize_Char_Arr(char* old_Char_Arr, long old_size, long new_size, bool is_First) {
        // create new array ptr
        char* new_Char_Arr = nullptr;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "entered resize char array function\n";
            cout << "New size: " << new_size << endl;
        }

        // attempt to allocate memory for new char array with updated size
        try{
            new_Char_Arr = new char[new_size];

        // catch any bad alloc errors and throw them up to file handler function
        } catch (const std::bad_alloc& e) {
            // display error message for current function and return failure
            cerr << "Memory allocation failed in resize_Char_Arr: " << e.what() << endl;
            //cout << "Memory allocation failed in resize_Char_Arr:\n";
            return nullptr;
        }

        long index;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "Memory resize for char array was sucessful, adding old data into new array\n";
        }

        // check and make sure this is not the first creation of the array
        if (!is_First && old_Char_Arr != nullptr) {

            // copy the old array into the new array index by index
            // for (index = 0; index < old_size; index++) {
            //     new_Char_Arr[index] = old_Char_Arr[index];
            // }

            //strcpy_s(new_Char_Arr, new_size, old_Char_Arr);
            copy_String(new_Char_Arr, new_size, genome_Size, old_Char_Arr );
    

            // toggle for displaying debug statements
            if( debug_Statements ) {
                cout << "Deleting old array\n";
            }

            // free the memory of the old array
            delete[] old_Char_Arr;
        }

        // return updated array
        return new_Char_Arr;
    }

    string* resize_Str_Arr(string* old_Str_Arr, long old_size, long new_size, bool is_First) {
        // create new array ptr
        string* new_Str_Arr = nullptr;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "entered resize str array function\n";
            cout << "new size: " << new_size << endl;
        }

        // attempt to allocate memory to new array of strings with updated size
        try{
            new_Str_Arr = new string[new_size];

        // capture any bad allocation errors and pass them back to file handler function
        } catch (const std::bad_alloc& e) {
            // display error message for current funciton and return failure
            cerr << "Memory allocation failed in resize_Str_Arr: " << e.what() << endl;
            //cout << "Memory allocation failed in resize_Str_Arr:\n";
            return nullptr;
        }

        long index;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "Memory resize for str array was sucessful, adding old data into new array\n";
        }

        // check and make sure this is not the first creation of the array
        if (!is_First && old_Str_Arr != nullptr) {
            // copy the contents from the old array over to the new array
            for (index = 0; index < old_size; index++) {

                //cout << "current name being copied from old str array: " << old_Str_Arr[index] << endl;
                new_Str_Arr[index] = old_Str_Arr[index];
            }

            // toggle for displaying debug statements
            if( debug_Statements ){
                cout << "Deleting old array\n";
            }

            // free the memory of the old array
            delete[] old_Str_Arr;
        }

        // return updated array
        return new_Str_Arr;
    }

    int* resize_Int_Arr(int* old_Int_Arr, long old_size, long new_size, bool is_First) {
        // create new array of Ints and set to nullptr
        int* new_Int_Arr = nullptr;
        
        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "entered resize int array function\n";
            cout << "new size: " << new_size << endl;
        }

        // attempt to allocate memory with updated size
        try{
            new_Int_Arr = new int[new_size];

        } catch (const std::bad_alloc& e) {
            cerr << "Memory allocation failed in resize_Int_Arr: " << e.what() << endl;
            //cout << "Memory allocation failed in resize_Int_Arr:\n";
            return nullptr;
        }

        int index;

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "Memory resize for int array was sucessful, adding old data into new array\n";
        }

        // check and make sure this is not the first creation of the array
        if (!is_First && old_Int_Arr != nullptr) {
            // copy the contents from the old array over to the new array
            for (int index = 0; index < old_size; index++) {
                new_Int_Arr[index] = old_Int_Arr[index];
            }

            // toggle for displaying debug statements
            if( debug_Statements ){
                cout << "Deleting old array\n";
            }

            // free the memory of the old array
            delete[] old_Int_Arr;
        }

        // return updated array
        return new_Int_Arr;
    }

    // merge sort function
    void mergeSort( char** arr, int low, int high){
        if( low < high ){
            int mid = low + ((high - low)/2);

            mergeSort(arr, low, mid);
            mergeSort(arr, mid+1, high);
            merge(arr, low, mid, high);
        }
    }

    // this will need to be able to sort the queries into alphabetical order
        // check the first letter
        // move the first letter either left or right depending on the alphabetical order
            // if the first letter is the same check the next letter
            // continue until the characters are different and can be sorted
    void merge(char** arr, int low, int mid, int high){
        // initialize variables
        int indexLeft, indexRight, kVal, leftLength, rightLength;
        leftLength = mid-low+1;
        rightLength = high - mid;
        char leftArr[leftLength][fragment_Size], rightArr[rightLength][fragment_Size];
        int char_Index;

        // add left side values to temp left array
        for(indexLeft=0; indexLeft < leftLength; indexLeft++){
            // leftArr[indexLeft] = *arr[low+indexLeft];
            copy_Col(leftArr[indexLeft], arr[low+indexLeft]);
        }


        // add right side values to temp right array
        for(indexRight=0; indexRight < rightLength; indexRight++){
            // rightArr[indexRight] = *arr[mid+1+indexRight];
            copy_Col(rightArr[indexRight], arr[mid+indexRight]);
        }

        // reset index values
        indexLeft = 0, indexRight = 0, kVal = 1;

        // sort the arrays !!!!!!!!!Need to add in logic for testing alphabetical position
        while( indexLeft < leftLength && indexRight < rightLength){
            // check if the value at the left index is less than the value on the right index
            if(leftArr[indexLeft] <= rightArr[indexRight]){
                // arr[kVal] = leftArr[indexLeft];
                copy_Col(arr[kVal], leftArr[indexLeft]);
                indexLeft++;
            }
            // otherwise assume right value is less
            else{
                //arr[kVal] = rightArr[indexRight];
                copy_Col(arr[kVal], rightArr[indexRight]);
                indexRight++;
            }
            kVal++;
        }

        // add the left side of the array to the return array
        while(indexLeft<leftLength){
            // arr[kVal] = leftArr[indexLeft];
            copy_Col(arr[kVal], leftArr[indexLeft]);
            kVal++;
            indexLeft++;
        }

        // add the right side of the array to the return array
        while(indexRight<rightLength){
            // arr[kVal] = rightArr[indexRight];
            copy_Col(arr[kVal], rightArr[indexRight]);
            kVal++;
            indexRight++;
        }
    }

    void bubbleSort(char** arr, int len){
        int innerIndex, outerIndex;
        bool swapped;

        for(outerIndex = 0; outerIndex < len-1; outerIndex++){
            swapped = false;
            for(innerIndex = 0; innerIndex < len-outerIndex-1; innerIndex++){
                if( arr[innerIndex] > arr[innerIndex+1]){
                    swap(arr[innerIndex], arr[innerIndex+1]);
                    swapped=true;
                }
            }

            if(!swapped){
                break;
            }
        }
        
    }

    void swap(char *oneVal, char *otherVal){

        int temp = *oneVal;
        *oneVal = *otherVal;
        *otherVal = temp;
    }

    long linear_Search(string target, char** query, long start, long end, long query_Length){
        int index; 
        
        for(index=0; index < query_Length; index++){
            if(target == query[index]){
                return index;
            }
        }

        return -1;
    }

    long binary_Search(string target, char** query, long low, long high){
        while(low <= high){
            long mid = low =((high-low)/2);

            if(query[mid] == target){
                return mid;
            }

            else if( query[mid] < target ){
                low = mid+1;
            }
            else{
                high = mid-1;
            }
        }
        
        return -1;
    }
    
    void copy_Col(char dest[], char fragment[]){
        int index;

        for(index=0;index<fragment_Size;index++){
            dest[index] = fragment[index];
        }
    }

}; 

int main(int argc, char* argv[]){
    Queries_AR my_Query;
    
    
    // toggle for displaying debug statements
    //if(my_Query.debug_Statements){
    cout << "Debug mode on, program start\n";
    //}

    my_Query.initial_Construct();
    
    cout << "inital_Cunstruct function returned\n";
    
    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 4 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }

    
    cout << "Calling to read_Query function\n";
    
    // read in the entire query dataset and store it in an instace of the Querie_AR class
    if(my_Query.read_Qurey(argv[2])){
        
        cout << "Sucessful return from read_Query function\n";
        
        // check for if the program will be donig binary search based off the command line parameters
        if(argv[3] == "-binary"){
        
            cout << "About to begin sorting queries\n";
        
            // sort the query in alphabetical order then
            my_Query.sortFragments(my_Query.query_Data, my_Query.query_Size, true);
            
            cout << "Query sort completed\n";
        }
        
        
        // read in the entire subject dataset into a single, concatenated character array
        // if(file_reader(argv[1])){
        //     // search the queries for every fragment ( pass in flag for desiered search algorithm )

        //     // display search time for the first 5k, 10k, 100k, and 1M 32 character long fragments of the subject dataset within the query dataset

        //     // how long would it take to search for every possible 32-long character fragment of the subject dataset within the query dataset

        //     // display the first 15 fragments of the subject dataset along with its indicies that you found within the query_AR object
        // }


        // display program end
        cout << "Displaying the first 15 SORTED queries" << endl;

        for(int index=0; index<15; index++){
            for(int col; col<my_Query.fragment_Size;col++){
                cout << my_Query.query_Data[index][col] << endl;

            }
        }
    }
    
    cout << "Program End\n";
    
    return 0;
}
