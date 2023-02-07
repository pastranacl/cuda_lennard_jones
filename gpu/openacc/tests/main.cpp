#include <iostream>


using namespace std;



int main()
{	
	int N=1e10;
	double sum=0;
	
	cout << "Init seque " << endl;
	for(int i=0; i<N; i++) {
		sum += i;
	}

	cout << "End seque. Start parallel " << endl;

	#pragma acc parallel reduction(+:sum)
	for(int i=0; i<N; i++) {
		sum += i;
	}
	return 0;
}
