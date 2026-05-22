#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdbool.h>

#define CHUNK_SIZE 16384

bool** adj_matrix;
int* index_link;
int graph_nodes = 4000;

// Allocating dynamic memory for 2D array
bool** getAllocate(int n) {
    int i;
    bool** arr = (bool**)malloc((n) * sizeof(bool*));
    for (i = 0; i < n; i++) {
        arr[i] = (bool*)malloc((n) * sizeof(bool));
    }
    return arr;
}

// Deallocating dynamic memory for 2D array
void getDeAllocate(int n, int** arr) {
    int i;
    for (i = 0; i < n; i++) {
        free(arr[i]);
    }
    free(arr);
}

void readAndDecompressGzipFile(const char *filename) {
    gzFile file;
    char buffer[CHUNK_SIZE];
    int bytesRead;
    int idx_counter = 0;
    int value1_present = 0;
    int value2_present = 0;
    int u, v;
    adj_matrix = getAllocate(graph_nodes);

    // Open the gzipped file
    file = gzopen(filename, "rb");
    // if (!file) {
    //     fprintf(stderr, "Error opening file %s!\n", filename);
    //     return 1;
    // }

    // Read and decompress the file chunk by chunk
    do {
        bytesRead = gzread(file, buffer, CHUNK_SIZE);
        if (bytesRead < 0) {
            fprintf(stderr, "Error reading from file %s!\n", filename);
            gzclose(file);
            // return 1;
        }

        // Process the buffer, split lines by '\n'
        char *token;
        char *rest = buffer;
        while ((token = strtok_r(rest, "\n", &rest))) {
            // Split each line by ',' to extract first two values
            char *value1_str = strtok(token, ",");
            char *value2_str = strtok(NULL, ",");

            if (value1_str && value2_str) {
                // Convert string values to integers
                int value1 = atoi(value1_str);
                int value2 = atoi(value2_str);
                value1_present = value2_present = 0;
                for(int id = idx_counter; id >= 0; id--){
                    if(index_link[id] == value1){
                        value1_present = 1;
                        u = id;
                    }
                    if(index_link[id] == value2){
                        value2_present = 1;
                        v = id;
                    }
                    if(value1_present && value2_present)
                        break;
                }

                if(!value1_present){
                    index_link[idx_counter] = value1;
                    u = idx_counter;
                    idx_counter++;
                }
                if(!value2_present){
                    index_link[idx_counter] = value2;
                    v = idx_counter;
                    idx_counter++;
                }
                adj_matrix[u][v] = 1;
                adj_matrix[u][u] = 1;
                adj_matrix[v][v] = 1;


                // Process value1 and value2 as integers
                printf("Value 1: %d, Value 2: %d\n", value1, value2);
                if (idx_counter > 10)
                    break;
            }
        }
    } while (bytesRead > 0);

    // Close the file
    gzclose(file);

    // return adj_matrix;
}

int main() {
    const char *filename = "datasets/soc-sign-bitcoinalpha.csv.gz";
    index_link = (int*)malloc(graph_nodes * sizeof(int));
    
    // if ((adj_matrix = readAndDecompressGzipFile(filename)) != 0) {
    //     fprintf(stderr, "Failed to read and decompress file %s\n", filename);     
    // }
    readAndDecompressGzipFile(filename);
    printf("Printing Adjacency Matrix. \n");
    for(int i = 0; i < 10; i++){
            for(int j = 0; j < 10; j++){
                printf("%d ", adj_matrix[i][j]);
            }
            printf("\n");
    }
    printf("Printing index_link. \n");
    for(int j = 0; j < graph_nodes; j++){
        printf("%d ", index_link[j]);
    }
    return 0;
}
