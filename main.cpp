// includes
#include "main.h"
using namespace std;

// Queries_AR class
class Queries_AR {
    public:
    // storage for query dataset (2D Array)
    char **query_Data;
    long query_Size;
    long allocated_Size;

    // storage for subject dataset
    char* genome_Data;
    long genome_Size;
    long allocated_Genome_Size;
    
    long scaffold_Count = 0;

    int* found_Frags;
    
    const bool debug_Statements = false;
    const int fragment_Size = 33;
    
    // 5k, 10K, 100K, and 1M time stamps
    time_t prog_Start = 0,
           five_Thousand = 0,
           one_Hund_Thousand = 0,
           one_Million = 0,
           prog_End = 0;
           
    //const long MAX_FRAGMENTS = 125000000; 
    const long MAX_FRAGMENTS = 125000000; 
    const long MAX_SCAFFOLD_COUNT = 608;

    // default constructor function
    void initial_Construct(){
        query_Data = nullptr;
        genome_Data = nullptr;

        query_Size = 0;
        allocated_Size = 0;
        genome_Size = 0;

        found_Frags = nullptr;
    }

    // custom constructor function for resizing the mulitdimesional array
    void query_Constructor(const string new_Query){
        // check if the query size is larger than the allocated space for the array
        if(query_Size >= allocated_Size){
            // increment the size by 1 to save space complexity
            long new_Size = query_Size + 1;
            
            // create new pointer array with new size
            char** new_Data = new char*[new_Size];
            
            // copy pointers from old array to new array
            for(long index = 0; index < query_Size; index++){
                new_Data[index] = query_Data[index];
            }
            
            // free old array memory
            delete[] query_Data;
            
            // set query_Data pointer to new array
            query_Data = new_Data;
            
            // set new size for tracking
            allocated_Size = new_Size;
        }
        
        // create a new memory storage for the new fragment
        query_Data[query_Size] = new char[fragment_Size];
        
        // copy the data into the 33 character space array
        strncpy(query_Data[query_Size], new_Query.c_str(), fragment_Size);
        
        // add terminating character to the end of the fragment
        query_Data[query_Size][fragment_Size-1] = '\0';
        
        // incremnt fragment count
        query_Size++;
    }

    // read query dataset file function
    bool read_Qurey(const string& file_Name){
        ifstream file(file_Name);
        string current_Line;
        int query_Num = 0;
        time_t stop_Watch = 0;
        
        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }
        
        // error handling for mainly bad_allocs or segmentation faults
        try {
            // allocate the inital array storage, to save time complexity set to max frags
            allocated_Size = MAX_FRAGMENTS;
            query_Data = new char*[allocated_Size];
            
            // iterate through file contents line by line
            while ( query_Size < MAX_FRAGMENTS && getline(file, current_Line)) {
                // check for stop watch times and set times
                if(query_Size == 1000 - 1 || query_Size == 10000 - 1 ||
                   query_Size == 500000 - 1 || query_Size == 1000000 - 1){
                    time(&stop_Watch);
                    cout << " Read in " << query_Size + 1 << " fragments at: " << 
                                                        ctime(&stop_Watch) << endl;
                }
                
                // check if we are at a valid fragment
                if(current_Line[0] != '>' && !current_Line.empty()){
                    // increment query count
                    query_Num++;
                    
                    
                    if(debug_Statements){
                          cout << "rd_Qry-Q_NUM: " << query_Num << endl;
                          cout << "rd-Qry-curLn: " << current_Line << endl;
                    } 
                    
                    // resize query array and add new fragment to next free index
                    query_Constructor(current_Line);
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

        // return failure up to main
        return false;
    }

    // search function
    void searchQuery(bool linear){
        // if search strategy is linear
        if(linear){
            linear_Search();
            
        }
        // Binary Search?
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
    // function for reading genome file ( LOOK AT THIS WHEN WE GET THE INPUT FILE )
    bool file_reader(const string& file_Name) {
        // variables
        ifstream file(file_Name);
        string current_Scaffold_Name;
        string current_Line;
        string temp_Genome = "";
        long prev_genome_size;
        unsigned long index;
        time_t stop_Watch = 0;

        // check if the file can be opened
        if (!file.is_open()) {
            // display error and exit
            cerr << "Error: Unable to open file " << file_Name << "\n";
            return false;
        }

        // error handling for mainly bad_allocs or segmentation faults
        try {
            // initialize starting array, set to 3 billion to save on runtime complexity
            allocated_Genome_Size = 3000000000;
            genome_Data = new char[allocated_Genome_Size];
            
            // iterate through file contents line by line
            while ( scaffold_Count < MAX_SCAFFOLD_COUNT && getline(file, current_Line)) {
                // check if we hit a header
                if(current_Line[0] == '>' && !temp_Genome.empty()){
                    
                    // increment scaffold count by 1
                    scaffold_Count++;
                    
                    // check and set timestamps
                    if(scaffold_Count == 100 || 
                       scaffold_Count == 300 || scaffold_Count == 500){
                           
                        time(&stop_Watch);
                        cout << "  Read in first " << scaffold_Count << " scaffolds at: " << ctime(&stop_Watch) << endl; 
                    
                    }
                    
                    // calculate new size for genome array
                    prev_genome_size = genome_Size;
    
    
                    // toggle for displaying debug statements
                    if( debug_Statements ){
                        cout << "Resizing genome array to fit current scaffold\n";
                    }
                
                    // resize the genome array to fit last scaffold and copy data into new array
                    genome_Constructor(temp_Genome, prev_genome_size);
    
                    // toggle for displaying debug statements
                    if( debug_Statements ){
                        cout << "Resize completed adding new values to genome\n";
                    }
    
                    // reset the temp genome string to empty
                    temp_Genome = "";
                }

                else if (!current_Line.empty() && current_Line[0] != '>') {
                    // Iterate through current line and append characters to temp_Genome
                    for (index = 0; index < current_Line.length() && index < 80; index++) {
                        if(  current_Line[index] == 'A' 
                          || current_Line[index] == 'C' 
                          || current_Line[index] == 'G' 
                          || current_Line[index] == 'T' 
                          || current_Line[index] == 'N'){
                      
                       
                            // Append character to temp_Genome
                            temp_Genome += current_Line[index];
                        }
                    }
                }
            }

            // post append information of last scaffold
            if (!temp_Genome.empty()) {
                
                // toggle for displaying debug statements
                if( debug_Statements ){
                    cout << "Last scaffold post loop being processed\n";
                }
                
                // calculate new genome array size
                prev_genome_size = genome_Size;
                
                // toggle for displaying debug statements
                if( debug_Statements ) {
                    cout << "Resizing genome arr for last scaffold\n";
                }
                
                // resize the genome array to fit last scaffold and copy data into new array
                genome_Constructor(temp_Genome, prev_genome_size);
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

        // iterate through array of pointers and free each index
        for(index=0;index<size;index++){
            delete[] arr_To_Destroy[index];
        }

        // free main array
        delete[] arr_To_Destroy;
        arr_To_Destroy=nullptr;
    }

    // deconstructor for deallocating the genome character array
    void genome_Deconstructor(char* genome_arr){
        delete[] genome_arr;
        genome_arr=nullptr;
    }

    // supporting function for copying from one string to another at a specific index
    void copy_String(char* main_Str, long location, long total_Size, const string to_Add){
        long index = 0;
        long location_Index;
        
        // copy data from to_Add inot main_Str starting at the location parameter
        for(location_Index = location; location_Index < total_Size; location_Index++){
            main_Str[location_Index] = to_Add[index];
            index++;
        }
    }

    void genome_Constructor(const string to_Add, long cat_Index) {
        long index;
        long length = to_Add.length();
        
        try{
            // check if the array needs to be resized
            if(genome_Size + length >= allocated_Genome_Size){
                // calculate the new size of the genome
                long new_Size = genome_Size + length;
                
                // allocate new array size
                char* new_Genome_Arr = new char[new_Size];
                
                // copy data from old array to new one
                for(index = 0; index < genome_Size; index++){
                    new_Genome_Arr[index] = genome_Data[index];
                }
                
                // deallocate memory from old array
                delete[] genome_Data;
                
                // set genome pointer to new array
                genome_Data = new_Genome_Arr;
                
                // set the new allocated size for genome array
                allocated_Genome_Size = new_Size;
            }
        }
        catch (const std::bad_alloc& e){
            // send back error message and return failure
            cerr << "Memory allocation failed: " << e.what() << endl;
        } 
        // catch any other errors that occur and throw them up to main funciton
        catch(const std::exception& e){
            // send back error message and return failure
            cerr << "Exception occurred: " << e.what() << endl;
        }
            
        for(index = 0; index < length; index++){
            //add in the new data 
            genome_Data[cat_Index + index] = to_Add[index];
        }
        
        // set the new genome size
        genome_Size = cat_Index + to_Add.length();
    }

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
            //cout << "Memory allocation failed in resize_Str_Arr:\n";
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
            //cout << "Memory allocation failed in resize_Int_Arr:\n";
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
        // check if the low is still less than the high
        if( low < high ){
            // calculate the mid
            int mid = low + (high - low)/2;

            // break down the left side recursivly
            mergeSort( low, mid);
            
            // break down the right side recursivly
            mergeSort( mid+1, high);
            
            // sort the two sides together alphabetically
            merge( low, mid, high);
        }
    }

    // this will need to be able to sort the queries into alphabetical order
        // check the first letter
        // move the first letter either left or right depending on the alphabetical order
            // if the first letter is the same check the next letter
            // continue until the characters are different and can be sorted
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
        
        // set index values
        indexLeft = 0;
        indexRight = 0;
        kVal = low; 
        
        // iterate until one of the arrays are empty
        while(indexLeft < leftLength && indexRight < rightLength){
            // check if the left string comes before the right in alphabeitcal order
            if(strcmp(leftArr[indexLeft], rightArr[indexRight]) <= 0){
                // copy over to return array and increment index
                copy_Col(leftArr[indexLeft], true, kVal);
                indexLeft++;
            } 
            
            // othwewise assume right came before
            else {
                // copy over to return array and incremnt index    
                copy_Col(rightArr[indexRight], true, kVal);
                indexRight++;
            }
            
            // incremnt return array index
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

    void bubbleSort(char** arr, int len){
        int innerIndex, outerIndex;
        bool swapped;

        // outer loop that will run for N length or if the array is sorted
        for(outerIndex = 0; outerIndex < len-1; outerIndex++){
            // set trigger to off
            swapped = false;
            
            // loop through and push bigest value to the rightmost index
            for(innerIndex = 0; innerIndex < len-outerIndex-1; innerIndex++){
                // if the current index comes after the right swap places
                if( arr[innerIndex] > arr[innerIndex+1]){
                    swap(arr[innerIndex], arr[innerIndex+1]);
                    
                    // set trigger to on
                    swapped=true;
                }
            }

            // if array was sorted before outer loop ends, break the loop early
            if(!swapped){
                break;
            }
        }
    }

    // swap pointers to to values
    void swap(char *oneVal, char *otherVal){
        int temp = *oneVal;
        *oneVal = *otherVal;
        *otherVal = temp;
    }

    void linear_Search(){
        long index = 0; 
        long queries_Found = 0;
        long high, query_Itr;
        int query_Index;
        string target = "";
        long target_Index;
        
        // iterate through the genome array unil end of genome or all fragments were found
        while(index < genome_Size && queries_Found < query_Size){
            
            // set high to the end of a test fragment
            high = index + fragment_Size - 1;
            
            // check if the end of the test fragment is within scope
            if( high < genome_Size ){
                // reset the target string
                target = "";

                // read test fragemnt into test string
                for(target_Index = index; target_Index<high; target_Index++){
                   target += genome_Data[target_Index];
                }
            }
            
            // iterate over the index + 31 characters
            for( query_Itr = 0; query_Itr < query_Size; query_Itr++ ){
                // compare with the fragment 
                // if a fragment was fond then the compare func will return its index in the query array
                query_Index = compare_Query(target, query_Itr);
                
                // if that index thats returned is not -1 than change that index to value 1
                if( query_Index != -1 ){
                    found_Frags[query_Index] = 1;
                    
                    // increment the queries_Found count
                    queries_Found++;
                    
                    // check and set timestamps
                    if(queries_Found == 5000){
                        time(&five_Thousand);
                    }
                    else if(queries_Found == 100000){
                        time(&one_Hund_Thousand);
                    }
                    else if(queries_Found == 1000000){
                        time(&one_Million);
                    }
                    
                    // stop searching 
                    break;
                }
                
                // if we hit the stop number for sort end early
                if( queries_Found >= 1000000 ){
                    break;
                }
            }
            // increment index
            index++;
        }
    }

    void binary_Search(){
        long top_Range, index = 0, char_Index = 0;
        long high = query_Size;
        long low = 0;
        long queries_Found = 0;
        bool frag_Found = false;
        long target_Index, queries_Checked = 0;
        char target_Arr[fragment_Size]; 
        
        // loop through genome until end of genome or until all fragments are found
        while(index < genome_Size && queries_Found < query_Size && queries_Found < 1000000){
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
                    //   target += genome_Data[target_Index];
                    target_Arr[char_Index] = genome_Data[target_Index];
                    char_Index++;
                }
                
                // set end of fragment
                target_Arr[fragment_Size-1] = '\0';
            }
            
            // now check for current test fragment in query array
            // loop until low is at or greater the hgih
            while(low <= high && !frag_Found && queries_Checked < query_Size-1){
                // calculate the mid index
                long mid = low +((high-low)/2);
    
                // check if the mid index fragment is a match to the test fragment
                // also make sure the mid index was not already found
                if(strcmp(query_Data[mid], target_Arr) == 0 && found_Frags[mid] != 1){
                    // set the mid to found
                    found_Frags[mid] = 1;
                    
                    // increment the numebr of queries found
                    queries_Found++;
                    
                    // check and set timestamps
                    if(queries_Found == 5000){
                        time(&five_Thousand);
                    }
                    else if(queries_Found == 100000){
                        time(&one_Hund_Thousand);
                    }
                    else if(queries_Found == 1000000){
                        time(&one_Million);
                    }
                    
                    // set trigger to true
                    frag_Found = true;
                }
    
                // otherwise check if the test fragment comes after mid index fragment
                else if( strcmp(query_Data[mid], target_Arr) < 0 ){
                    // set the low to the right side min
                    low = mid+1;
                }
                // otherwise assume fragment comes before
                else{
                    // set high the max left side
                    high = mid-1;
                }
                
                // increment queries checked counter
                queries_Checked++;
            }
            index++;
        }
    }
    
    void copy_Col( char fragment[], bool to_Query, int row){
        int index;

        // copy data to the query 
        if(to_Query){
            for(index=0;index<fragment_Size;index++){
                query_Data[row][index] = fragment[index];
            }
            
        }
        
        // copy data from the query
        else{
            for(index=0;index<fragment_Size;index++){
                fragment[index] = query_Data[row][index];
            }
        }
    }
    
    int compare_Query(string target, long query_Row){
        int index;
        
        // iterate over the test fragment
        for( index = 0; index < fragment_Size-1; index++){
            
            // check if there is a differnce
            if(target[index] != query_Data[query_Row][index]){
                // return failure if difference
                //return query_Row;
                return -1;
            }
        }
        
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
    
    time(&stop_Watch);
    cout << "Program Start at: " << ctime(&stop_Watch) << endl;
    
    // toggle for displaying debug statements
    if(my_Query.debug_Statements){
        cout << "Debug mode on, program start\n";
        }

    // zero out the query and genome data
    my_Query.initial_Construct();
    
    // check if the program has the appropriate number of parameter
    // create error handler for no input file argument
    if( argc < 4 ){
        cout << "Error please input the correct command in\n Program End";
        return 1;
    }
    
    time(&stop_Watch);
    cout << "Beginning to read query file at: " << ctime(&stop_Watch) << endl;
    
    // read in the entire query dataset and store it in an instace of the Querie_AR class
    if(my_Query.read_Qurey(argv[2])){
        
        cout << "Sucessful return from read_Query function\n";
        cout << "There were " << my_Query.query_Size << " Fragments" << endl;
        
        cout << "\nDisplaying the first 15 unsorted queries" << endl;

        for(int index=0; index<15; index++){
            
            cout << " ";
            
            for(int col=0; col<my_Query.fragment_Size;col++){
                cout << my_Query.query_Data[index][col];

            }
            cout << endl;
        }
        
        // check for if the program will be donig binary search based off the command line parameters
        if(strcmp(argv[3], "-binary") == 0){
        
            time(&stop_Watch);
            cout << "Beginning fragment sort at: " << ctime(&stop_Watch) << endl;
        
            // turn binary trigger on
            binary_Search = true;
            
            // sort the query in alphabetical order then
            my_Query.sortFragments(my_Query.query_Data, my_Query.query_Size-1, true);
            
            // display program end
            cout << "\nDisplaying the first 15 SORTED queries" << endl;
    
            for( index=0; index<15; index++){
                
                cout << " ";
                
                for( col=0; col<my_Query.fragment_Size;col++){
                    cout << my_Query.query_Data[index][col];
    
                }
                cout << endl;
            }
        }
        
        time(&stop_Watch);
        cout << "\n Starting genome file reader at: " << ctime(&stop_Watch) << endl;
        
        
        // read in the entire subject dataset into a single, concatenated character array
        if(my_Query.file_reader(argv[1])){
            
            cout << "\nGenome size: " << my_Query.genome_Size << endl;

            // create fragment array of ints with each index set to 0
            my_Query.found_Frags = my_Query.resize_Int_Arr(my_Query.found_Frags, 0, 
             my_Query.query_Size, true);
            
            time(&stop_Watch);
            cout << "Starting Search at: " << ctime(&stop_Watch) << endl;
            
            // search the queries for every fragment ( pass in flag for desiered search algorithm )
            if( binary_Search ){
                cout << "Conducting binary search\n\n";
                my_Query.searchQuery(false);
            }
            else {
                cout << "Conducting Linear search\n\n";
                my_Query.searchQuery(true);
            }
            
            time(&stop_Watch);
            cout << "Finished Search at: " << ctime(&stop_Watch) << endl;
            // display search time for the first 5k, 10k, 100k, and 1M 32 character long fragments of the subject dataset within the query dataset

            // how long would it take to search for every possible 32-long character fragment of the subject dataset within the query dataset

            // display the first 15 fragments of the subject dataset along with its indicies that you found within the query_AR object
        
            cout << "Displaying first 15 fragments\n";
            for(index=0; index<15; index++){
                cout << " " << my_Query.found_Frags[index] << endl;
            }
            
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
    
    // deallocate the memory for both arrays
    my_Query.qurey_Deconstructor(my_Query.query_Data, my_Query.query_Size);
    my_Query.genome_Deconstructor(my_Query.genome_Data);
    
    cout << "\nProgram End\n";
    
    return 0;
}