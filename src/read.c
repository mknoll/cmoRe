#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef int elem_type ;

//http://ndevilla.free.fr/median/median/src/quickselect.c
#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type quick_select(elem_type arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

void median(int *n, int *x, int *res) {
    res[0] = quick_select(x, n[0]);
    //printf("%f\m", binapprox(6, val));
    //printf("%f\n", res[0]);
}

void sum(int *n, int *x, int *res) {
    res[0] = 0;
    for (int i=0; i < n[0]; i++) {
	res[0] = res[0]+x[i];
    }
    printf("RESULT: %d\n", res[0]);
    //printf("%f\m", binapprox(6, val));
    //printf("%f\n", res[0]);
}

//const char* getfield(char* line, int num)
const char* getfield(char* line, int num, int del)
{
    char *tok;    
    //char delimiter[1] = "\t";    
    char d1[1]= "\t";
    char d2[1] = ",";
    char d3[1] = ";";
    char *delimiter;
    if (del == 0) {
	delimiter = d2;
    } else if( del == 1) {
	delimiter = d1;
    } else if( del == 2) {
	delimiter = d3;
    }

    tok = strsep(&line, delimiter); // first token    
    while(tok != NULL) {        
	//printf ("%s ", tok);    
	// do something    
	if (!--num)        
	    return tok;    
	tok = strsep(&line, delimiter); /* next token */    
    }    
    return NULL;
}

// Get lines
void getNRows(char **file, int *n) {
    printf("opening: %s\n", file[0]);

    // get number of lines
    FILE* stream = fopen(file[0], "r");
    int lines = 0;
    while (EOF != (fscanf(stream, "%*[^\n]"), fscanf(stream,"%*c")))
	++lines;
    printf("Lines : %i\n", lines);
    fclose(stream);

    n[0] = lines;
}

void readCol(char **file, int *nrow, int *ncol, double *ret, int *del) {
    printf("opening: %s\n", file[0]);

    //get data
    FILE* stream = fopen(file[0], "r");
    char line[32768]; //TODO: check
    int i = 0;
    while (fgets(line, 32768, stream) != NULL)
    {
        char* tmp = strdup(line);
	//printf("Field 3 would be %s\n", getfield(tmp, 3));
        // NOTE strtok clobbers tmp
	if (i >= 1) {
	    ret[i] = atof(getfield(tmp, ncol[0], del[0]));
	    //printf("%f\n", ret[i]);
	} 
	i++;
        free(tmp);
    }
    fclose(stream);
}


void readColChar(char **file, int *nrow, int *ncol, char **ret, int *del) {
    printf("opening: %s\n", file[0]);

    //get data
    FILE* stream = fopen(file[0], "r");
    char line[32768]; //TODO: check
    int i = 0;
    while (fgets(line, 32768, stream) != NULL)
    {
        char* tmp = strdup(line);
	//printf("Field 3 would be %s\n", getfield(tmp, 3));
        // NOTE strtok clobbers tmp
	if (i >= 1) {
	    const char *stringA = getfield(tmp, ncol[0], del[0]);
	    strncpy(ret[i], stringA, strlen(stringA)+1);
	    //printf("%f\n", ret[i]);
	} 
	i++;
        free(tmp);
    }
    fclose(stream);
}
