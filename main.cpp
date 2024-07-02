// includes
#include "main.h"
using namespace std;

class Queries_AR {
    public:
    // storage for query dataset (2D Array)
    char **query_Data;
    long query_Size;
    int* found_Frags;

    // storage for subject dataset
    char* genome_Data;
    long genome_Size;
    long scaffold_Count = 0;
    
    const bool debug_Statements = false;
    
    // 5k, 10K, 100K, and 1M time stamps
    time_t prog_Start = 0,
           five_Thousand = 0,
           one_Hund_Thousand = 0,
           one_Million = 0,
           prog_End = 0;
           
    const int fragment_Size = 33;
    const long MAX_FRAGMENTS = 100000000; 
    const long MAX_SCAFFOLD_COUNT = 10000;

    // default constructor function
    void initial_Construct(){
        query_Data = nullptr;
        genome_Data = nullptr;

        query_Size = 0;
        genome_Size = 0;

        found_Frags = nullptr;
    }

    // custom constructor function for resizing the mulitdimesional array
    void query_Constructor(const int new_Size, const bool is_First, const string new_Query){
        int index, row, col;
        char **new_Temp = nullptr;

        // allocate initial single array of pointers
        new_Temp = new char *[new_Size];
    
        // loop thorugh initial array and allocate spaces for each index needed
        for(index=query_Size;index<new_Size;index++){
            // give the pointer some memory
            new_Temp[index] = new char[fragment_Size];
        }

        // check if this is NOT the first query
        if(!is_First){
            // copy data into new array ( at same index )
            for(row=0;row<query_Size;row++){
                new_Temp[row] = query_Data[row];
            }
            
            // delete old array
            delete[] query_Data;
        }

        // add new data into the new spot
        for(index=0;index<fragment_Size-1;index++){
            new_Temp[new_Size-1][index] = new_Query[index];
        }
        
        // add null terminating character
        new_Temp[new_Size-1][fragment_Size-1] = '\0';

        // set pointer to new array
        query_Data = new_Temp;
    }

    // read query dataset file function
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string current_Line;
        int query_Num = 0;
        time_t stop_Watch = 0;
        bool runLoop = true;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        // error handling for mainly bad_allocs or segmentation faults
        try {
            // iterate through file contents line by line
            while ( query_Size < MAX_FRAGMENTS && getline(file, current_Line)) {
                
                // timestamps for checkng program progress
                if(query_Size == 1000 - 1){
                    time(&stop_Watch);
                    cout << " Read in 1000 fragments at: " << ctime(&stop_Watch) << endl;
                }
                
                if(query_Size == 10000 - 1){
                    time(&stop_Watch);
                    cout << " Read in 10000 fragments at: " << ctime(&stop_Watch) << endl;
                }
                
                if(query_Size == 500000 - 1){
                    time(&stop_Watch);
                    cout << " Read in 500000 fragments at:" << ctime(&stop_Watch) << endl;
                }
                
                if(query_Size == 1000000 - 1){
                    time(&stop_Watch);
                    cout << " Read in 1000000 fragments at: " << ctime(&stop_Watch) << endl;
                }
                
                // check if we are at a valid fragment
                if(current_Line[0] != '>' && !current_Line.empty()){
                    // increment query count
                    query_Num++;
                    
                    if(debug_Statements){
                          cout << "rd_Qry-Q_NUM: " << query_Num << endl;
                          cout << "rd-Qry-curLn: " << current_Line << endl;
                    } 

                    // resize query array and add new data
                    query_Constructor(query_Num, query_Num == 1, current_Line);

                    // increment query size variable in class
                    query_Size++;
                }

                // check if the end of file has been reached to kill the loop
                if(file.eof()){
                    cout << "End of file reached, breaking loop\n";
                    break;
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

        // return failure
        return false;
    }

    // search function
    void searchQuery(bool linear){
        // if linear search was chosen
        // O(n)
        if(linear){
            linear_Search();    
        }
        // otherwise assume binary Search
        // O(query_Size log query_Size)
        else{
            binary_Search();
        }
    }

    // sort fragments function
    void sortFragments(char** query_Data, int size, bool ms){
        // sort all character fragments in alpahbetical order (lexicographic) order (logorithmic runtime!)
        // check for merge sort code
        if(ms){
            mergeSort( 0, size);
        }
        // otherwise assume bubble sort
        else{
            bubbleSort(query_Data, size);
        }
    }

    // function to read in the genome data
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string current_Scaffold_Name;
        string current_Line;
        string temp_Genome = "";
        long new_genome_size;
        long protein_Count = 0;
        long genome_Index;
        unsigned long index;
        time_t stop_Watch = 0;
        bool runLoop = true;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }

        // error handling for mainly bad_allocs or segmentation faults
        try {
            // iterate through file contents line by line
            while ( scaffold_Count < MAX_SCAFFOLD_COUNT && getline(file, current_Line)) {
                
                // check if we hit a header
                if(current_Line[0] == '>' && !temp_Genome.empty()){
                    
                    // checks for program progress
                    if(scaffold_Count == 100){
                        time(&stop_Watch);
                        cout << "  Read in first 100 scaffolds at: " << ctime(&stop_Watch) << endl; 
                    }
                    
                    if(scaffold_Count == 300){
                        time(&stop_Watch);
                        cout << "  Read in first 300 scaffolds at: " << ctime(&stop_Watch) << endl; 
                    }
                    
                    if(scaffold_Count == 500){
                        time(&stop_Watch);
                        cout << "  Read in first 500 scaffolds at: " << ctime(&stop_Watch) << endl; 
                    }
                    
                    // increment scaffold count by 1
                    scaffold_Count++;
                    
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
                    genome_Data = resize_Char_Arr(genome_Data, new_genome_size, 
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

                // otherwise assume valid character stream for genome
                else if (!current_Line.empty() && current_Line[0] != '>') {
                    // Iterate through current line and append characters to temp_Genome
                    for (index = 0; index < current_Line.length() && index < 80; index++) {
                        if( current_Line[index] == 'A' 
                         || current_Line[index] == 'C' 
                         || current_Line[index] == 'G' 
                         || current_Line[index] == 'T'
                         || current_Line[index] == 'N'){

                            // Append character to temp_Genome
                            temp_Genome += current_Line[index];
                        }
                    }
                }
                
                if(file.eof()){
                    cout << "End of file reached, breaking loop\n";
                    break;
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

                // resize the genome array to fit last scaffold and copy data into new array
                genome_Data = resize_Char_Arr(genome_Data, new_genome_size, genome_Size == 0);

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

        // free memory for each pointer inside the array
        for(index=0;index<size;index++){
            delete[] arr_To_Destroy[index];
        }

        // free memory for main array and set to null ptr
        delete[] arr_To_Destroy;
        arr_To_Destroy=nullptr;
    }


    // deconstructory for genome character array
    void genome_Deconstructor(char* genome_arr){
        // free the memory for the character array and set pointer to nullptr
        delete[] genome_arr;
        genome_arr=nullptr;
    }


    // copy to the contents of a string over to a character array
    void copy_String(char* main_Str, long location, long total_Size, const string to_Add){
        long index = 0;
        long location_Index;
        
        // loop through the character array starting at a specified position and add the string characters
        for(location_Index = location; location_Index < total_Size; location_Index++){
            main_Str[location_Index] = to_Add[index];
            index++;
        }
    }


    // function for resizing the genome character array
    char* resize_Char_Arr(char* old_Char_Arr, long new_size, bool is_First) {
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
            return nullptr;
        }

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "Memory resize for char array was sucessful, adding old data into new array\n";
        }

        // check and make sure this is not the first creation of the array
        if (!is_First && old_Char_Arr != nullptr) {
            // copy the data from the old array over to the new one
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


    // function for resizing the array mainly responsible for scaffolds 
    string* resize_Str_Arr(string* old_Str_Arr, long old_size, long new_size, bool is_First) {
        // create new array ptr
        string* new_Str_Arr = nullptr;
        long index;

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
            return nullptr;
        }

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


    // function for resizing an array of integers
    int* resize_Int_Arr(int* old_Int_Arr, long old_size, long new_size, bool is_First) {
        // create new array of Ints and set to nullptr
        int* new_Int_Arr = nullptr;
        int index;
        
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
            return nullptr;
        }

        // toggle for displaying debug statements
        if( debug_Statements ){
            cout << "Memory resize for int array was sucessful, adding old data into new array\n";
        }

        // check and make sure this is not the first creation of the array
        if (!is_First && old_Int_Arr != nullptr) {
            // copy the contents from the old array over to the new array
            for (index = 0; index < old_size; index++) {
                new_Int_Arr[index] = old_Int_Arr[index];
            }

            // toggle for displaying debug statements
            if( debug_Statements ){
                cout << "Deleting old array\n";
            }

            // free the memory of the old array
            delete[] old_Int_Arr;
        }

        else{
            // zero out indexs
            for (index = 0; index < new_size; index++) {
                    new_Int_Arr[index] = 0;
                }
        }
        
        // return updated array
        return new_Int_Arr;
    }

    // merge sort function
    void mergeSort( int low, int high){
        if( low < high ){
            int mid = low + (high - low)/2;

            // sort the left side first recursively
            mergeSort( low, mid);

            // sort right side next recursivley
            mergeSort( mid+1, high);

            // merge sort left and right side together
            merge( low, mid, high);
        }
    }

    
    // merge left array and right array together in a alphabetically sorted order
    void merge( int low, int mid, int high){
        // initialize variables
        int indexLeft, indexRight, kVal, leftLength, rightLength;
        leftLength = mid-low+1;
        rightLength = high - mid;
        char leftArr[leftLength][fragment_Size], rightArr[rightLength][fragment_Size];

        // Copy left side values to temp left array
        for(indexLeft = 0; indexLeft < leftLength; indexLeft++){
            copy_Col(leftArr[indexLeft], false, low + indexLeft);
        }

        // Copy right side values to temp right array
        for(indexRight = 0; indexRight < rightLength; indexRight++){           
            copy_Col(rightArr[indexRight], false, mid + 1 + indexRight);
        }

        // set indexs for arrays
        indexLeft = 0;
        indexRight = 0;
        kVal = low; 

        // loop thorugh until one array has no more items left
        while(indexLeft < leftLength && indexRight < rightLength){

            // check if left character is smaller than right character
            if(strcmp(leftArr[indexLeft], rightArr[indexRight]) <= 0){
                // copy left character into main array and increment left index
                copy_Col(leftArr[indexLeft], true, kVal);
                indexLeft++;
            } 
            
            // otherwise assume right was smaller
            else {  
                // copy right character into main array and increment right index      
                copy_Col(rightArr[indexRight], true, kVal);
                indexRight++;
            }

            // increment main array index
            kVal++;
        }
        
        // Copy remaining elements from leftArr
        while(indexLeft < leftLength){
            copy_Col(leftArr[indexLeft], true, kVal);
            indexLeft++;
            kVal++;
        }
        
        // Copy remaining elements from rightArr
        while(indexRight < rightLength){ 
            copy_Col(rightArr[indexRight], true, kVal);
            indexRight++;
            kVal++;
        }
    }


    // sort array by pushing largest item to right most side of the array
    void bubbleSort(char** arr, int len){
        int innerIndex, outerIndex;
        bool swapped;

        // main loop for sort
        for(outerIndex = 0; outerIndex < len-1; outerIndex++){
            // set trigger to off
            swapped = false;

            // loop thorugh and push large items to the right
            for(innerIndex = 0; innerIndex < len-outerIndex-1; innerIndex++){
                // if left character is larger than right
                if( arr[innerIndex] > arr[innerIndex+1]){
                    // swap indexs
                    swap(arr[innerIndex], arr[innerIndex+1]);

                    // set trigger to on
                    swapped=true;
                }
            }

            // if no swap occured in inner loop then sort is finished
            if(!swapped){
                break;
            }
        }
    }


    // swap the pointers of two items
    void swap(char *oneVal, char *otherVal){
        int temp = *oneVal;
        *oneVal = *otherVal;
        *otherVal = temp;
    }

    
    // search for an item in an array starting from left side all the way to the right
    void linear_Search(){
        long index = 0; 
        long queries_Found = 0;
        long high, query_Itr;
        long query_Index;
        string target = "";
        long target_Index;
        
        // iterate through the genome array unil end of genome or all fragments were found
        while(index < genome_Size && queries_Found < query_Size){
            // calculate end of fragment in genome
            high = index + fragment_Size - 1;
            
            // check if its within the genome array
            if( high < genome_Size ){
                // reset target string
                target = "";

                // copy fragment into target test string
                for(target_Index = index; target_Index<high; target_Index++){
                   target += genome_Data[target_Index];
                }
            }
            
            // check each query starting from starting index to ending index
            for( query_Itr = 0; query_Itr < query_Size; query_Itr++ ){
                // compare with the fragment 
                // if a fragment was fond then the compare func will return its index in the query array
                query_Index = compare_Query(target, query_Itr);
                
                // if that index thats returned is not -1 than change that index to value 1
                if( query_Index != -1 ){
                    found_Frags[query_Index] = 1;
                    
                    // increment the queries_Found count
                    queries_Found++;

                    // set timestamps
                    if(queries_Found == 5000){
                        time(&five_Thousand);
                    }
                    else if(queries_Found == 100000){
                        time(&one_Hund_Thousand);
                    }
                    else if(queries_Found == 1000000){
                        time(&one_Million);
                    }
                }
            }  
            
            index++;         
        }
    }


    // search for a fragment in divide and conquerer approach
    void binary_Search(){
        long top_Range, index = 0, char_Index = 0;
        long high = query_Size;
        long low = 0;
        long queries_Found = 0;
        bool frag_Found = false;
        long target_Index, queries_Checked = 0;
        char target_Arr[fragment_Size]; 
        
        // loop through genome until end of genome or until all fragments are found
        while(index < genome_Size && queries_Found < query_Size){
            // reset variables for new loop
            frag_Found = false;
            queries_Checked = 0;
            low = 0;
            high = query_Size-1;
            
            // create the high for the test fragment
            top_Range = index + fragment_Size - 1;
            
            // test if the high is in bounds
            if( top_Range < genome_Size ){
                char_Index = 0;

                // read in 32 characters from genome ( index to top_Range )
                for(target_Index = index; target_Index<top_Range; target_Index++){
                    target_Arr[char_Index] = genome_Data[target_Index];
                    char_Index++;
                }
                
                // add terminating character to end to match query data
                target_Arr[fragment_Size-1] = '\0';
            }
            
            // now check for current test fragment in query array
            // loop until low is at or greater the hgih
            while(low <= high && !frag_Found && queries_Checked < query_Size-1){
                // calculate the middle
                long mid = low +((high-low)/2);
    
                // check if the mid value was the fragment and has not already been found
                if(strcmp(query_Data[mid], target_Arr) == 0 && found_Frags[mid] != 1){
                    
                    // set fragment to found
                    found_Frags[mid] = 1;

                    // increment found fragment counter
                    queries_Found++;
                    
                    // set timestamps
                    if(queries_Found == 5000){
                        time(&five_Thousand);
                    }
                    else if(queries_Found == 100000){
                        time(&one_Hund_Thousand);
                    }
                    else if(queries_Found == 1000000){
                        time(&one_Million);
                    }
                    
                    // set trigger to on
                    frag_Found = true;
                }
    
                // othwerwise check if the mid is less than the target (alphabetical order wise)
                else if( strcmp(query_Data[mid], target_Arr) < 0 ){
                    // set low to right side 
                    low = mid+1;
                }

                // otherwise assume on left side
                else{
                    // set high to left side
                    high = mid-1;
                }
                
                // increment checked counter to avoid over checking
                queries_Checked++;
            }
            index++;
        }
    }
    

    // copy data to and from the query array
    void copy_Col( char fragment[], bool to_Query, int row){
        int index;

        // check if we need to copy data into query array
        if(to_Query){
            for(index=0;index<fragment_Size;index++){
                query_Data[row][index] = fragment[index];
            }
        }

        // otherwise assume taking data out of query array
        else{
            for(index=0;index<fragment_Size;index++){
                fragment[index] = query_Data[row][index];
            }
        }
    }
    

    // compare if the query and the fragment are the same
    int compare_Query(string target, long query_Row){
        int index;
        
        // loop through the fragment string
        for( index = 0; index < fragment_Size-1; index++){

            // if the character does not match
            if(target[index] != query_Data[query_Row][index]){              
                //return query_Row;
                return -1;
            }
        } 

        // otherwise assume match      
        //return -1;
        return query_Row;
    }
};

int main(int argc, char* argv[]){
    Queries_AR my_Query;
    int index, col;
    bool binary_Search = false;
    time(&my_Query.prog_Start);
    time_t stop_Watch = 0;
    
    // stamp program start time
    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;
    
    // toggle for displaying debug statements
    if(my_Query.debug_Statements){
        cout << "Debug mode on, program start\n";
    }

    // initialize query data
    my_Query.initial_Construct();
    
    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 4 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }
    
    // timestamp the start of the query file reader
    time(&stop_Watch);
    cout << "Beginning to read query file at: " << ctime(&stop_Watch) << endl;
    
    // read in the entire query dataset and store it in an instace of the Querie_AR class
    if(my_Query.read_Qurey(argv[2])){
        
        // program progress outputs
        cout << "Sucessful return from read_Query function\n";
        cout << "There were " << my_Query.query_Size << " Fragments" << endl;
        cout << "\nDisplaying the first 15 unsorted queries" << endl;

        // loop through first 15 queries and display to output
        for(int index=0; index<15; index++){
            cout << " ";

            for(int col=0; col<my_Query.fragment_Size;col++){
                cout << my_Query.query_Data[index][col];
            }

            cout << endl;
        }
        
        // check for if the program will be donig binary search based off the command line parameters
        if(strcmp(argv[3], "-binary") == 0){
        
            // timestamp for sort start
            time(&stop_Watch);
            cout << "Beginning fragment sort at: " << ctime(&stop_Watch) << endl;
        
            // if binary flag was true then set trigger for binary sort
            binary_Search = true;

            // sort the query in alphabetical order then
            my_Query.sortFragments(my_Query.query_Data, my_Query.query_Size-1, true);
            
            // display first 15 sorted queries
            cout << "\nDisplaying the first 15 SORTED queries" << endl;
            for( index=0; index<15; index++){
                cout << " ";
                
                for( col=0; col<my_Query.fragment_Size;col++){
                    cout << my_Query.query_Data[index][col];
    
                }

                cout << endl;
            }
        }
        
        // timestamp for genome file reader
        time(&stop_Watch);
        cout << "\n Starting genome file reader at: " << ctime(&stop_Watch) << endl;
        
        // read in the entire subject dataset into a single, concatenated character array
        if(my_Query.file_reader(argv[1])){
            
            cout << "\nGenome size: " << my_Query.genome_Size << endl;

            // create fragment array of ints with each index set to 0
            my_Query.found_Frags = my_Query.resize_Int_Arr(my_Query.found_Frags, 0, 
                                                          my_Query.query_Size, true);
            
            // timestamp for search start
            time(&stop_Watch);
            cout << "Starting Search at: " << ctime(&stop_Watch) << endl;
            
            // search the queries for every fragment ( pass in flag for desiered search algorithm )
            if( binary_Search ){
                cout << "Conducting binary search\n\n";
                my_Query.searchQuery(false);
            }

            // otherwise assume linear search strategy
            else {
                cout << "Conducting Linear search\n\n";
                my_Query.searchQuery(true);
            }
            
            // timestamp for search completion
            time(&stop_Watch);
            cout << "Finished Search at: " << ctime(&stop_Watch) << endl;
            
            // display the first 15 fragments of the subject dataset along with its indicies that you found within the query_AR object
            cout << "Displaying first 15 fragments\n";
            for(index=0; index<15; index++){
                cout << " " << my_Query.found_Frags[index] << endl;
            }
            
            // display all timestamps
            cout << "\n\n";
            cout << "PROGRAM STARTED     : " << ctime(&my_Query.prog_Start) << endl;
            
            if( my_Query.five_Thousand != 0 ){
                cout << "First 5k fragments  : " << ctime(&my_Query.five_Thousand) << endl;
            }
            if( my_Query.one_Hund_Thousand != 0 ){
                cout << "First 100k fragments: " << ctime(&my_Query.one_Hund_Thousand) << endl;
            }
            if( my_Query.one_Million != 0 ){
                cout << "First 1M fragments  : " << ctime(&my_Query.one_Million) << endl;
            }
            
            time(&my_Query.prog_End);
            cout << "PROGRAM END         : " << ctime(&my_Query.prog_End) << endl;
            
        }
    }
    
    // free memory for query and genome arrays
    my_Query.qurey_Deconstructor(my_Query.query_Data, my_Query.query_Size);
    my_Query.genome_Deconstructor(my_Query.genome_Data);
    
    cout << "\nProgram End\n";
    
    return 0;
}