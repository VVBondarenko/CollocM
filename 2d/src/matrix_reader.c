#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
	FILE *matr;
	matr =fopen("matrix","r");
	int i,j,N = atoi(argv[1]);
	float cache;
	for(i = 0; i<N; i++)
	{
		for(j = 0; j<N; j++)
		{
			fscanf(matr,"%f", &cache);
			printf("%3.2g\t", cache);
		}
		printf("\n");
	}

		
	return 0;
}
