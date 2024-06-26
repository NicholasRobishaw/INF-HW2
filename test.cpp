#include "main.h"

using namespace std;

int main(){
    string str_One = "Big";
    string str_Two = "Apples";
    string str_Three = "Bid";

    string return_Arr[3] = {str_One, str_Two, str_Three};

    cout << "Initial array unsorted\n";

    for(int index=0; index < 3; index++){
        cout << return_Arr[index] << endl;
    }

    bubbleSort(return_Arr, 3);

    cout << "\nInitial array sorted\n";
    
    for(int index=0; index < 3; index++){
        cout << return_Arr[index] << endl;
    }

}

void bubbleSort(string arr[], int len){
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

void swap(int *oneVal, int *otherVal){
    int temp = *oneVal;
    *oneVal = *otherVal;
    *otherVal = temp;
}