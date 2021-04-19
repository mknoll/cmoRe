#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

const char* getfield(char* line, int num)
{
    char *tok;    
    char delimiter[1] = "\t";    
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

void readCol(char **file, int *nrow, int *ncol, double *ret) {
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
	    ret[i] = atof(getfield(tmp, ncol[0]));
	    //printf("%f\n", ret[i]);
	} 
	i++;
        free(tmp);
    }
    fclose(stream);
}


void readColChar(char **file, int *nrow, int *ncol, char **ret) {
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
	    const char *stringA = getfield(tmp, ncol[0]);
	    strncpy(ret[i], stringA, strlen(stringA)+1);
	    //printf("%f\n", ret[i]);
	} 
	i++;
        free(tmp);
    }
    fclose(stream);
}
