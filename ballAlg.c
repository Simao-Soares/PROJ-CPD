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
    long id;
    struct _node *L;
    struct _node *R;
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
node_t *build_tree(double **pts, int n_dims,int n_points);
void dump_tree(node_t *root);
double **computeMostDistantPoints(double **pts, int n_dims, int n_points);
double computeDistance(double *p1, double *p2, int n_dims, int starter);
void orthogonalProjection(double **pts, double *a,double *b,int n_dims, int n_points);
int findMedian(double **pts, int n_points,int n_dims, double *a);


/* main: process parameters */
int main(int argc, char *argv[])
{
  double exec_time;
  exec_time = -omp_get_wtime();
  double **pts;
  int n_dims = 0;
  long n_points = 0;
  node_t *root;

  pts = get_points(argc, argv, &n_dims, &n_points);
  root = build_tree(pts, n_dims, n_points);
  exec_time += omp_get_wtime();
  fprintf(stderr, "%.1lf\n", exec_time);
  dump_tree(root); // to the stdout!
}


void dump_tree(node_t *root)
{

}

//calculate the distance between two points
double computeDistance(double *p1, double *p2, int n_dims, int starter)
{
    double result = 0;
    int i = 0;

    for(i=0; i<n_dims; i++)
    {
      printf("%f - %f\n", p1[i], p2[i+starter]);
      result = result + ((p1[i] - p2[i+starter]) * (p1[i] - p2[i+starter]));
    }
    result = sqrt(result);
    printf("%f\n\n", result);
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


}

//choose the two points that have a greater distance
double **computeMostDistantPoints(double **pts, int n_dims, int n_points)
{
  int i = 0;
  int a_index = 0 ;
  double greaterDistance = 0;
  double possibleGreaterDistance = 0;
  double **mostDistant = (double **)malloc(2 * sizeof(double*));
  for(i = 0; i < 2; i++)
  {
    mostDistant[i] = (double *)malloc(n_dims * sizeof(double));
  }

  for(i=1; i<n_points ; i++){
    possibleGreaterDistance = computeDistance(pts[0], pts[i], n_dims, 0);
    if( possibleGreaterDistance > greaterDistance)
      {
        greaterDistance = possibleGreaterDistance;
        a_index = i;
        mostDistant[0] = pts[i];
      }
  }
  for(i=0; i<n_points && i != a_index; i++){ // ver se a condiçao do i afeta a performance
    possibleGreaterDistance = computeDistance(pts[a_index], pts[i], n_dims, 0);
    if( possibleGreaterDistance > greaterDistance)
      {
        greaterDistance = possibleGreaterDistance;
        mostDistant[1] = pts[i];
      }  
  }

  /*for(i=0; i < n_points; i++)
  {
    for(j = i+1; j< n_points; j++)
    {
      possibleGreaterDistance = computeDistance(pts[i], pts[j], n_dims, 0);
      if( possibleGreaterDistance > greaterDistance)
      {
        greaterDistance = possibleGreaterDistance;
        mostDistant[0] = pts[i];
        mostDistant[1] = pts[j];
      }
    }
  }*/
  return mostDistant;
}


int findMedian(double **pts, int n_points,int n_dims, double *a)
{
  int i = 0;
  int j =0;
  int sumLefts = 0;
  int medianIndex = n_points/2;
  int notLeaf = 1;
  qNode *root = NULL;
  root = (qNode *) malloc(sizeof(qNode));
  root->right = NULL;
  root->left = NULL;
  root->nLefts = 0;
  qNode *aux = NULL;
  aux = (qNode *) malloc(sizeof(qNode));

  for ( i = 0; i < n_points; i++)
  {
    for ( j = 0; j < n_dims; j++)
    {
      printf("ponto %d : %f \n",i,pts[i][j]);
    }
    printf("\n");
    
  }
  

  
  root->value = computeDistance(a, pts[medianIndex], n_dims, n_dims);
  root->id = medianIndex;

  for(i=0; i< n_points; i++)
  {
    if(i == medianIndex)
    {
      continue;
    }
    for(j=0;j<n_dims;j++)
    {
      printf("%d!-%f\n",i, pts[i][j]);
    }
    printf("\n");
    qNode *me = NULL;
    me = (qNode *) malloc(sizeof(qNode));
    me->right = NULL;
    me->left = NULL;
    me->value = computeDistance(a, pts[i], n_dims, n_dims);
    me->nLefts = 0;
    me->id = i;
    aux = root;
    notLeaf = 1;
    while(notLeaf)
    {
      if(me->value > aux->value){
        if(aux->right == NULL){
          aux->right = me;
          notLeaf = 0;
        }
        else{
          aux = aux->right;
        }
      }
      else{
        aux->nLefts ++;
        if(aux->left == NULL){
          aux->left = me;
          notLeaf = 0;
        }
        else{
          aux = aux->left;
        }
      }
    }
  }

  aux = root;
  while(aux->nLefts + sumLefts != medianIndex)
  {

    if(aux->nLefts + sumLefts < medianIndex){
      if(aux->right != NULL){
        sumLefts = sumLefts + aux->nLefts + 1;
        aux = aux->right;
      }
      else{
        return aux->id;
      }
    }
    else{
      if(aux->left != NULL){
        aux = aux->left;
      }
      else{
        return aux->id;
      }
    }
  }
  return aux->id;

}


//create the ball tree
node_t *build_tree(double **pts, int n_dims,int n_points)
{
  int i=0;
  int index = 0;
  double **mostDistant = (double **)malloc(2 * sizeof(double*)); //two most distant points in a set (matrix)
  for(i = 0; i < 2; i++)
  {
    mostDistant[i] = (double *)malloc(n_dims * sizeof(double));
  }

  mostDistant = computeMostDistantPoints(pts, n_dims, n_points);
  computeDistance(mostDistant[0], mostDistant[1], n_dims, 0);
  orthogonalProjection(pts,mostDistant[0],mostDistant[1], n_dims, n_points);

  /*for ( i = 0; i < n_points; i++)
  {
    for ( int j = 0; j < n_dims; j++)
    {
      printf("ponto %d : %f \n",i,pts[i][j]);
    }
    printf("\n");
  }*/
  
  index = findMedian(pts, n_points, n_dims, mostDistant[0]);
  printf("median:\n");
  for(i=0; i < n_dims;i++)
  {
    printf("%f\n",pts[index][i+n_dims]);
  }




  return NULL;
}

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * 2 *sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * 2 * n_dims];

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
        for(j = 0; j < *n_dims; j++)
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;

#ifdef DEBUG
    for(i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}
