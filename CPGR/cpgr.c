/*Implementation of cpgr*/

/* 

Use:  gcc cpgr.c -lm -o cpgr -O1 -lrt -lz for compilation. Change the optimization level as required.

*/

/*On windows follow this instruction to get a running time:
  1.clock_t start = clock();
  2.clock_t stop = clock();
  3.elapsed = ((double)(stop - start)) / CLOCKS_PER_SEC * 1000.0;
 */

 /*On Grid follow this instruction to get a running time:
   0. struct timespec begin, end;
   1.clock_gettime(CLOCK_REALTIME, &begin);
   2.clock_gettime(CLOCK_REALTIME, &end);         
   3.long seconds = end.tv_sec - begin.tv_sec;
   4.long nanoseconds = end.tv_nsec - begin.tv_nsec;
   5.elapsed = (seconds + nanoseconds * 1e-9) * 1000;
  */

/*

for Segmentation Fault Error Please Check the following variable paths:
1: f_name
2: saveFile 
3: tempFile in Function save_graph_to_mtx()


*/

#pragma warning(disable:4996)
#define _CRT_SECURE_NO_WARNINGS

#include<stdio.h>
#include<string.h>
#include<stdbool.h>
#include<stdint.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <zlib.h>

#define MAXCHAR 100000
#define CHUNK_SIZE 16384
#define split true
FILE* k_values;   // Stores k_hat values
FILE* m_hat_values;
FILE* results;
FILE* saveFile;
FILE* tempFile;

struct timespec begin, end;
char f_name[150]; // name of the adjacency matrix file
float compression_ratio;
int** adj_matrix;
int* d_v;  // Degree of vertices in W
int* K;  // Stores initial selection of vertices from W
int* K_split;     // Stores the right partition for each delta-clique
int* U_split;      // Stores the left partition for each delta-clique
int clique_u_size;     // Size of the left partition of the delta-clique
int clique_v_size;      // Size of the right partition of the delta-clique
int k_split;    // starting index of q in each cpgr iteration
int k_temp;     // ending index of q in each cpgr iteration
int m_hat; // number of remaining edges
int d_K;  // degree of set K in iteration k
int Gamma;     // # of cliques extracted in each iteration of cpgr
int k_hat;
float delta;
int density;
int nodes;

char cores[] = "cpgr";
int multiplier;
char saveFilename[150];
int experiment_no;
int* temp_psi;  // stores the degrees of vertices w for sorting 
int* temp_psi_idx;  // stores the sorted indices of degrees of vertices w 
int initial_edges;
int total_edges;    // Sum fo edges in delta-clique and trivial edges
int graph_nodes;
double execution_time;
int* rightPartitionSize;
int* leftPartitionSize;
int middlePartitionSize = 0;
int* tempMiddlePartition;
int* edges;
int cliqueIndex; // Starting index of clique vertex
int edges_in_clique = 0;
bool saveToFile = true;



// Allocating dynamic memory for 2D array
int** getAllocate(int n) {
    int i;
    int** arr = (int**)malloc((n) * sizeof(int*));
    for (i = 0; i < n; i++) {
        arr[i] = (int*)malloc((n) * sizeof(int));
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



/* Algo2: cpgr lines 5 and 13 */
void get_k_hat() {
    float de = (2 * pow((double)graph_nodes, 2)) / m_hat;
    float nu = delta * log2((double)graph_nodes);
    k_hat = floor((double)nu / (log2((double)de)));
}

// To check the format of the given file (csv or mtx)
int csvFile() {
    int fileNameLen = strlen(f_name);
    int targetLen = 3;
    char fileType[] = "csv";
    if (fileNameLen < targetLen) {
        return 0; // Target is longer than the file name, not possible to match
    }

    for (int i = fileNameLen - targetLen, j = 0; i < fileNameLen; i++, j++) {
        if (f_name[i] != fileType[j]) {
            return 0; // Mismatch found, target not found in the last three characters of the file name
        }
    }

    return 1; // Target found in the last three characters of the file name
}

//To read the input from csv file
void load_adj_matrix() {
    adj_matrix = getAllocate(graph_nodes);
    FILE* fpointer = fopen(f_name, "r");
    char line[MAXCHAR];
    if (!fpointer)
        printf("Cann't open the file\n");
    else {
        int row = 0;
        int col = 0;
        while (fgets(line, MAXCHAR, fpointer)) {
            col = 0;
            //strtok break down each line into smaller string
            char* value = strtok(line, ",");
            while (value != NULL) {
                adj_matrix[row][col] = atoi(value);
                //adj_matrix_fix[row][col] = atoi(value);
                value = strtok(NULL, ",");
                col++;
            }
            row++;
        }
        fclose(fpointer);
    }
}


//to read the input from .mtx file
void readMatrixMarketFile() {
    FILE* file = fopen(f_name, "r");
    // Check if the file could be opened
    if (file == NULL) {
        printf("Failed to open the file.\n");
        exit(1);
    }

    // Read the header information
    char line[256];
    fgets(line, sizeof(line), file);
    if (strncmp(line, "%%MatrixMarket matrix coordinate", 31) != 0) {
        printf("Invalid Matrix Market file.\n");
        exit(1);
    }

    fgets(line, sizeof(line), file);
    while (line[0] == '%') {
        fgets(line, sizeof(line), file);
    }

    // Parse the graph size and number of edges
    sscanf(line, "%d %d %d \n", &leftPartitionSize, &rightPartitionSize, &edges);
    int i;
    // Setting graph nodes as a maximum nodes among left and right partition
    // We are also initializing the clique index  
    if (leftPartitionSize > rightPartitionSize) {
        cliqueIndex = leftPartitionSize;
        graph_nodes = leftPartitionSize;
    }
    else {
        cliqueIndex = rightPartitionSize;
        graph_nodes = rightPartitionSize;
    }
    // Initializing the dynamic array for adjacency matrix
    adj_matrix = getAllocate(graph_nodes);
    for (i = 0; i < edges; i++) {
        int row, col;
        fscanf(file, "%d %d \n", &row, &col); 
        adj_matrix[row - 1][col - 1] = 1;
    }
    m_hat = edges;
    fclose(file);
}

// reading .gz file
int readAndDecompressGzipFile() {
    gzFile file;
    char buffer[CHUNK_SIZE];
    int bytesRead;

    // Open the gzipped file
    file = gzopen(f_name, "rb");
    if (!file) {
        fprintf(stderr, "Error opening file %s!\n", f_name);
        return 1;
    }

    // Read and decompress the file chunk by chunk
    do {
        bytesRead = gzread(file, buffer, CHUNK_SIZE);
        if (bytesRead < 0) {
            fprintf(stderr, "Error reading from file %s!\n", f_name);
            gzclose(file);
            return 1;
        }
        // Process the buffer, e.g., write it to another file
        // Here we just print it to the console
        fwrite(buffer, 1, bytesRead, stdout);
    } while (bytesRead > 0);

    // Close the file
    gzclose(file);

    return 0;
}




// Temporarily saving edges of the extracted cliques
void saveCliquesEdges() {
    cliqueIndex++;
    for (int i = 0; i < clique_u_size; i++) {
        fprintf(tempFile, "%d %d\n", U_split[i], cliqueIndex);
    }
    for (int j = 0; j < clique_v_size; j++) {
        fprintf(tempFile, "%d %d\n", cliqueIndex, K_split[j]);
    }
}

// Saving the compressed graph in .mtx format
void save_graph_to_mtx() {
    // Write the header
    fprintf(saveFile, "%%MatrixMarket matrix coordinate pattern general\n");
    fprintf(saveFile, "%% Compressed graph edges\n");
    fprintf(saveFile, "%% -------------------------------------------\n");
    fprintf(saveFile, "%% Original Graph: %s\n", f_name);
    fprintf(saveFile, "%% Graph_nodes:%d, compression_ratio:%f\n", graph_nodes, compression_ratio);
    fprintf(saveFile, "%% -------------------------------------------\n");
    fprintf(saveFile, "%% leftPartitionSize, middlePartitionSize, rightPartitionSize, edges\n");
    fprintf(saveFile, "%d %d %d %d\n", leftPartitionSize, middlePartitionSize, rightPartitionSize, total_edges);
    
    for (int i = 0; i < leftPartitionSize; i++) {
        for (int j = 0; j < rightPartitionSize; j++) {
            if (adj_matrix[i][j])
                fprintf(saveFile, "%d %d\n", i + 1, j + 1);
        }
    }

    // Update the directory path as required
    // tempFile = fopen("datasets/tempCliqueEdges.mtx", "r");
    tempFile = fopen("datasets/tempCliqueEdges.mtx", "r");
    if (tempFile == NULL) {
        printf("Failed to open the file 265 .\n");
        exit(1);
    }

    // Read the header information
    char line[256];
    char ch;
    
    while (true) {
        int u, v;
        fscanf(tempFile, "%d %d \n", &u, &v); 
        if (u == -1 || v == -1)
            break;
        fprintf(saveFile, "%d %d\n", u, v);
    }
    fclose(tempFile);
}


/* Algo2: cpgr lines 1 and 11 */
void get_edges() {
    m_hat = 0;
    int i;
    for (i = 0; i < graph_nodes; i++)
        m_hat += d_v[i];
}

/* Used by partition function */
void swap(int* xp, int* yp) {
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}



// Partition the array and return the pivot index
int partition(int* arr, int* temp_psi_idx, int low, int high) {
    int pivot = arr[low]; // Choose the first element as the pivot
    int i = low + 1; // Start from the second element
    int j = high;

    while (i <= j) {
        // Find an element greater than the pivot from the left side
        while (i <= j && arr[i] >= pivot) {
            i++;
        }
        // Find an element smaller than the pivot from the right side
        while (i <= j && arr[j] < pivot) {
            j--;
        }
        // Swap the elements if they are out of order
        if (i < j) {
            swap(&arr[i], &arr[j]);
            swap(&temp_psi_idx[i], &temp_psi_idx[j]);
        }
    }

    // Move the pivot to its correct position
    swap(&arr[low], &arr[j]);
    swap(&temp_psi_idx[low], &temp_psi_idx[j]);

    return j;
}


/* Algo3: CSA line 2*/
// Quicksort implementation
void quicksort(int* arr, int* temp_psi_idx, int low, int high) {
    if (low < high) {
        int pivotIndex = partition(arr, temp_psi_idx, low, high);

        //#pragma omptask shared(arr, temp_psi_idx) firstprivate(low, pivotIndex)
        quicksort(arr, temp_psi_idx, low, pivotIndex - 1);

        //#pragma omptask shared(arr, temp_psi_idx) firstprivate(high, pivotIndex)	
        quicksort(arr, temp_psi_idx, pivotIndex + 1, high);
        //#pragma omptaskwait
    }
}


/* Algo3: CSA lines 4 to 6*/
void get_new_K(int k) {
    int i = 0;
    for (i = 0; i < graph_nodes; i++)
        temp_psi_idx[i] = i;
    memcpy(temp_psi, d_v, graph_nodes * sizeof(int));
    (void)quicksort(temp_psi, temp_psi_idx, 0, graph_nodes - 1);

    i = 0;
    //printf(" Right Partition \n");
    while (d_v[temp_psi_idx[i]] >= temp_psi[k_hat - 1]) {
        K[i] = temp_psi_idx[i] + 1;
        //printf("%d ", K[k][i]);
        i++;

    }
    K[i] = -1; // to define end of selection of vertex v in V
    d_K = i;
}

/* Algo3: CSA lines 12 and 13 */
void get_updates(int d) {
    //int d;
    int j, i;
    for (j = 0; j < clique_v_size; j++) {
        d_v[K_split[j] - 1] = d_v[K_split[j] - 1] - clique_u_size;
        m_hat -= clique_u_size;
        for (i = 0; i < clique_u_size; i++) {
            adj_matrix[U_split[i] - 1][K_split[j] - 1] = 0;
        }
    }
    edges_in_clique += clique_u_size + clique_v_size;
}



/* Algo3: CSA lines 8 to 11 */
void get_U_with_k_hat(int k) {

    int div_K = floor(((double)d_K / k_hat));
    int start = 0;
    int end = k_hat;
    int div = 0;
    int set_size = k_hat;
    k_temp = k_temp + div_K;
    int i, j, d;

    /* Algo3: CSA line 9 */

    for (d = k_split; d < k_temp; d++) {
        for (j = 0; j < set_size; j++) {
            K_split[j] = K[j + start];
        }
        K_split[set_size] = -1;
        clique_v_size = j;
        start = start + k_hat;
        
        /* Algo3: CSA line 10 */
        int u = 0;
        for (i = 0; i < graph_nodes; i++) {
            int flag = 0;
            for (j = 0; j < k_hat + 1 && !flag; j++) {
                if (K_split[j] != -1) {
                    if (adj_matrix[i][K_split[j] - 1] == 0) {
                        flag = 1;
                    }
                }
                else {
                    U_split[u] = i + 1;
                    u++;
                    flag = 1;
                }
            }
            U_split[u] = -1;
        }
        clique_u_size = u;
        get_updates(d);
        if (saveToFile) {
            saveCliquesEdges();
            middlePartitionSize++;
        }
    }
}

/* Algo2: cpgr lines 2 to 4*/
void get_d_v() {
    int j, i;
    for (j = 0; j < graph_nodes; j++) {
        d_v[j] = 0;
        for (i = 0; i < graph_nodes; i++) {
            if (adj_matrix[i][j] == 1)
                d_v[j] = d_v[j] + adj_matrix[i][j];
            else
                adj_matrix[i][j] = 0;
        }
    }
}


void display_adj_matrix() {
    printf("\nADJACENCY MATRIX\n");
    int i, j;
    for (i = 0; i < graph_nodes; i++) {
        for (j = 0; j < graph_nodes; j++) {
            printf(" %d  ", adj_matrix[i][j]);
        }
        printf("\n");
    }
}



/* Calculating compression ratio */
void get_compression_ratio() {
    get_edges();
    int edges_in_clique = 0;
    int p;
    for (p = 0; p < k_temp; p++) {
        if (clique_v_size < 2)
            edges_in_clique += clique_u_size;
        else
            edges_in_clique += clique_u_size + clique_v_size;
    }
    total_edges = m_hat + edges_in_clique;
    compression_ratio = (float)initial_edges / (float)total_edges;
}




void sequentialCPA() {
    clock_gettime(CLOCK_REALTIME, &begin);
    //clock_t start = clock();
    int k = 0;
    get_d_v();
    initial_edges = m_hat;
    get_k_hat();
    int flag = 0;
    while (k_hat > 1 && !flag)
    {
        int m = m_hat;
        get_new_K(k);

        get_U_with_k_hat(k);
        get_k_hat();

        k++;
        k_split = k_temp;
        if (m_hat == m) {
            flag = 1;
        }
    }



    //clock_t stop = clock();
    //run_time = ((double)(stop - start)) / CLOCKS_PER_SEC * 1000.0;
    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    execution_time = (seconds + nanoseconds * 1e-9) * 1000;
    // get_compression_ratio();
}



/* Run on windows without batch script */
/* updated main function to read mtx and csv files */
int main(int argc, char* argv[]) {
    if (argc < 5) {
        printf("Insufficient command-line arguments!\n");
        return 1;
    }

    nodes = atoi(argv[1]);
    density = atoi(argv[2]);
    experiment_no = atoi(argv[3]);
    delta =  atof(argv[4]);
    // Update the directory path as required
    // sprintf(f_name, "/ocean/projects/cis230093p/srabin/Graph_Compression/dataset1/bipartite_graph_%d_%d_%d.mtx", nodes, density, experiment_no);
    sprintf(f_name, "datasets/bipartite_graph_%d_%d_%d.mtx", nodes, density, experiment_no);
    const char *extension = strrchr(f_name, '.');
    
    
    if (strcmp(extension, ".gz") == 0) {
        return readAndDecompressGzipFile();
    } else if (csvFile()) {
        load_adj_matrix();
    } else {
        readMatrixMarketFile();
    }



    k_temp = 0;
    k_split = 0;
    // Update the directory path as required
    sprintf(saveFilename, "datasets/cpgr_tripartite_graph_%d_%d_%d_%d.mtx", nodes, density, experiment_no, (int) delta);
    // sprintf(saveFilename, "/ocean/projects/cis230093p/achavan/Graph_Compression/dataset/cpgr_tripartite_graph_%d_%d_%d_%d.mtx", nodes, density, experiment_no, (int)(delta*100));
    multiplier = ceil(log10((double)graph_nodes));
    d_v = (int*)malloc(graph_nodes * sizeof(int));
    K = (int*)malloc((graph_nodes + 1) * sizeof(int));
    K_split = (int*)malloc((graph_nodes + 1) * sizeof(int));
    U_split = (int*)malloc((graph_nodes + 1) * sizeof(int));
    // printf("%d, %d, %d, %f, %f, %f\n", graph_nodes, density, experiment_no, delta, compression_ratio, execution_time);
    saveFile = fopen(saveFilename, "w");
    // Update the directory path as required
    // tempFile = fopen("/ocean/projects/cis230093p/achavan/Graph_Compression/dataset/tempCliqueEdges.mtx", "w");
    tempFile = fopen("datasets/tempCliqueEdges.mtx", "w");
    temp_psi = (int*)malloc(graph_nodes * sizeof(int));
    temp_psi_idx = (int*)malloc(graph_nodes * sizeof(int));
    sequentialCPA();
    fprintf(tempFile, "%d", -1);
    fclose(tempFile);
    total_edges = m_hat + edges_in_clique;
    compression_ratio = (float)initial_edges / (float)total_edges;
    printf("%d, %d, %d, %f, %f, %f\n", graph_nodes, density, experiment_no, delta, compression_ratio, execution_time);
    save_graph_to_mtx();
    free(K);

    getDeAllocate(graph_nodes, adj_matrix);
    free(temp_psi);
    free(temp_psi_idx);
    fclose(saveFile);

    free(d_v);
    // Update the directory path as required
    // remove("/ocean/projects/cis230093p/achavan/Graph_Compression/dataset/tempCliqueEdges.mtx");
    remove("datasets/tempCliqueEdges.mtx");
    return 0;
}
