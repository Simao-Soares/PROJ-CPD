//Program that resolves The Ball Algorithm serially

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include<math.h>

#define RANGE 10
//#define DEBUG 1
#define random rand
#define srandom srand

typedef struct _node {
    double radius;
    int id;
    struct _node *L;
    struct _node *R;
    double *center;
} node_t;

//node of Q-Sort algorithm
typedef struct _qnode {
    double value;
    int id;
    int nLefts; //number of nodes under its left child
    struct _qnode *left;
    struct _qnode *right;
} qNode;

double **create_array_pts(int n_dims, long np);
double **get_points(int argc, char *argv[], int *n_dims, long *np);
node_t *build_tree(double **pts, int n_dims,int n_points, int *id);
void dump_tree(node_t *root,int n_dims);
double **computeMostDistantPoints(double **pts, int n_dims, int n_points);
double computeDistance(double *p1, double *p2, int n_dims, int starter);
void orthogonalProjection(double **pts, double *a,double *b,int n_dims, int n_points);
int findMedian(double **pts, int n_points,int n_dims, double *a);
void quickSort(double **array, int low, int high, int n_dims);
void swap(double **a, double **b);
int partition(double **array, int low, int high, int n_dims);


/* main: process parameters */
int main(int argc, char *argv[])
{
  double exec_time;
  exec_time = -omp_get_wtime();
  double **pts;
  double **pts_aux;
  double *_pts_aux;
  int n_dims = 0;
  long n_points = 0;
  node_t *root;
  int id = 0;

  pts = get_points(argc, argv, &n_dims, &n_points);
  pts_aux = pts ;
  _pts_aux = pts[0];
  #pragma omp parallel
  {
    #pragma omp single
      root = build_tree(pts_aux, n_dims, n_points, &id);
  }
  exec_time += omp_get_wtime();
  fprintf(stderr, "%.1lf\n", exec_time);
  printf("%d %d \n", n_dims, id+1);
  dump_tree(root,n_dims); // to the stdout!

  free(_pts_aux);
  free(pts);
}


void dump_tree(node_t *node,int n_dims)
{
  int i =0;
  
  if(node->L != NULL){
    printf("%d %d %d %f ",node->id, (node->L)->id,(node->R)->id, node->radius);
    for(i=0;i<n_dims;i++){
      printf("%f ",(node->center)[i]);
    }
    printf("\n");
    dump_tree(node->L, n_dims);
  }
  else{
    printf("%d -1 -1 %f ",node->id, node->radius);
    for(i=0;i<n_dims;i++){
      printf("%f ",(node->center)[i]);
    }
    printf("\n");
    //free(node->center);
    free(node);
    return;
  } 
  if(node->R != NULL){
    dump_tree(node->R, n_dims);
  }
  else{
    //free(node->center);
    //free(node);
    return;
  }
  free(node->center);
  free(node); 
}

//calculate the distance between two points
double computeDistance(double *p1, double *p2, int n_dims, int starter)
{
    double result = 0;
    int i = 0;
    
    for(i=0; i<n_dims; i++)
    {
      result = result + ((p1[i] - p2[i+starter]) * (p1[i] - p2[i+starter]));
    }
    result = sqrt(result);
    return result;
}

//perform the orthogonal projection of all points onto line ab
void orthogonalProjection(double **pts, double *a,double *b,int n_dims, int n_points)
{
  int i = 0;
  int j =0;
  double result1 = 0;
  double result2 = 0;
  double *b_a = (double *)malloc(n_dims * sizeof(double*));

  for(j=0;j<n_points;j++)
  {
    result1 = 0;
    result2 = 0;
    for(i = 0; i<n_dims; i++)
    {
      pts[j][n_dims + i]= pts[j][i]-a[i];
      b_a[i] = b[i] - a[i];
      result1 = result1 + (pts[j][n_dims +i] * b_a[i]);
      result2 = result2 + (b_a[i] * b_a[i]);
    }

    for(i = 0; i<n_dims; i++)
    {
      pts[j][n_dims + i] = ((result1/result2) * b_a[i]) + a[i];
    }
  }
  free(b_a);
}

//choose the two points that have a greater distance
double **computeMostDistantPoints(double **pts, int n_dims, int n_points)
{
  int i = 0;
  int initial_index = 0;
  int a_index = 0;
  double greaterDistance = 0;
  double possibleGreaterDistance = 0;

/* Allocation of A and B points */
  double **mostDistant = (double **)malloc(2 * sizeof(double*));
 

// Find the most Left Point 
  for(i=1;i<n_points;i++)
  {
    if(pts[initial_index][n_dims*2 +1] > pts[i][n_dims*2 + 1])
    {
      initial_index = i;
    }
  }
  
// Find Point A 
  for(i=0; i<n_points ; i++){
    possibleGreaterDistance = computeDistance(pts[initial_index], pts[i], n_dims, 0);
    if( possibleGreaterDistance > greaterDistance)
      {
        greaterDistance = possibleGreaterDistance;
        a_index = i;
        mostDistant[0] = pts[i];
      }
  }
  greaterDistance = 0;

// Find Point B
  for(i=0; i<n_points; i++){ 
    possibleGreaterDistance = computeDistance(pts[a_index], pts[i], n_dims, 0);
    if( possibleGreaterDistance > greaterDistance)
      {
        greaterDistance = possibleGreaterDistance;
        mostDistant[1] = pts[i];
      }  
  }
  return mostDistant;
}



void swap(double **a, double **b) {
  double *t = NULL;
  t = *a;
  *a = *b;
  *b = t;
}

// Function to partition the array on the basis of pivot element
int partition(double **array, int low, int high, int n_dims) {
  
  // Select the pivot element
  double pivotDistance = array[high][n_dims];
  int i = (low - 1);

  // Put the elements smaller than pivot on the left 
  // and greater than pivot on the right of pivot
  for (int j = low; j < high; j++) {
    if (array[j][n_dims] <= pivotDistance) {
      i++;
      swap(&array[i], &array[j]);
    }
  }
  swap(&array[i + 1], &array[high]);
  return (i + 1);
}
/* Quick Sort Algorithm */
void quickSort(double **array, int low, int high, int n_dims) {
  if (low < high) {
  
  // Select pivot position and put all the elements smaller 
  // than pivot on left and greater than pivot on right
  int pi = partition(array, low, high, n_dims);
  
  // Sort the elements on the left of pivot
  quickSort(array, low, pi - 1, n_dims);
  
  // Sort the elements on the right of pivot
  quickSort(array, pi + 1, high, n_dims);
  }
}


/************************************* BALL CREATION ******************************************/
node_t *build_tree(double **pts, int n_dims,int n_points,int *id)
{
  /* Ball Allocation */
  node_t *node = (node_t *)malloc(sizeof(node_t));

  /****** Leaf Case Creation*******/
  if(n_points == 1)
  {
    /* Radius, ID, L, R, Center Filling */
    node->radius = 0;
    node->id = *id ;
    node->L = NULL;
    node->R = NULL;
    node->center = pts[0];
    return node;
  }

  node->id = *id;


  /****** Not Leaf Creation *****/
  int i=0;
  //int j=0;
  int index = 0;
  double radius1 = 0;
  double radius2 = 0;
  int odd = n_points % 2;
  double **mostDistant = NULL; //two most distant points in a set (matrix)

  /* Find A and B points */
  mostDistant = computeMostDistantPoints(pts, n_dims, n_points);

  /* Orthogonal projetion points Creatrion */
  orthogonalProjection(pts,mostDistant[0],mostDistant[1], n_dims, n_points);

  //index = findMedian(pts, n_points, n_dims, mostDistant[0]);

  /* Sorting Points by distance */
  quickSort(pts,0,n_points-1, n_dims);
  index = n_points/2;

  /* node->center Allocation */
  double* center = (double *)malloc(n_dims * sizeof(double));

  /* Ball Center  Value and Assignment */
  for(i=0; i < n_dims;i++)
  {
    center[i] = (pts[index][i+ n_dims]+pts[index-1+odd][i+n_dims])/2 ;
  }

  node->center = center;

  /* Ball Radius computation  */
  if((radius1 = computeDistance(mostDistant[0],center,n_dims,0))>(radius2=computeDistance(mostDistant[0],center,n_dims,0))){
    node->radius = radius1;
  }
  else{
    node->radius = radius2;
  }

  free(mostDistant); 

  /* Node ID association */
      /* Left Ball Creation */
      #pragma omp task shared(id)
        {
          *id = *id +1; 
          node->L = build_tree(pts, n_dims, n_points/2, &(*id));
        }
      
      #pragma omp taskwait

      /* Right Ball Creation*/
      #pragma omp task shared(id)
        {
          *id = *id +1; 
          node->R = build_tree(&pts[n_points/2], n_dims, n_points/2 + odd, &(*id));
        }
  return node;
}


double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;
    long i = 0;

    _p_arr = (double *) malloc((n_dims * np * 2 + np) *sizeof(double));/* Allocation of points and space for orthogonal points */
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }
    //#pragma omp parallel for private(i)
    for(i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * 2 * n_dims + i];

    return p_arr;
}


double **get_points(int argc, char *argv[], int *n_dims, long *np)
{
    double **pt_arr;
    unsigned seed;
    long i;
    int j;

    if(argc != 4){
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    *n_dims = atoi(argv[1]);
    if(*n_dims < 2){
        printf("Illegal number of dimensions (%d), must be above 1.\n", *n_dims);
        exit(2);
    }

    *np = atol(argv[2]);
    if(*np < 1){
        printf("Illegal number of points (%ld), must be above 0.\n", *np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **) create_array_pts(*n_dims, *np);

    for(i = 0; i < *np; i++)
    {
        pt_arr[i][(*n_dims * 2) + 1] = i;
        for(j = 0; j < *n_dims; j++)
        {
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;
        }
    }

#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}