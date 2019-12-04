#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* globals ... */

#define HEADER "**********************************************************\n" \
               "             JACOBI ITERATION SOLUTION METHOD             \n" \
               "    (for strictly diagonally dominant square matrices)    \n" \
               "**********************************************************\n"
#define EPSILON 0.0000001
#define MAX_PRINT_ITERS 333

char SHOW_ITERS = 0; /* defaults to NO ... */

/* input matrix elements */
void getInput( int, double** mat );
/* check if diagonally dominant */
int checkIfDD( int, double** mat );
/* calculate and display answers */
void jacobiCalcDisplay( int, double** mat );

double dotProd( int num, double* v1, double* v2 );

void showXcheck( int num, double** mat, double* v1, double* v2 );

int checkFlags( int num, int flags[] );

/* free dynamic memory occupied by 2D square matrix */
void freeMat( int num, double** sqMat );

int more()/* defaults to yes ... */
{
    int reply ;
    printf("\nMore ... (y/n) ? ");
    reply = getchar();
    if( reply != '\n' )
        while( getchar() != '\n' ) ; /* flush stdin ... */
    putchar( '\n' );
    return !(reply=='n' || reply=='N');
}



int main()
{
    double **mat;
    int num;

    do
    {
        puts( HEADER );

        printf( "Show each iteration result 1=Yes, 0=No ? " );
        scanf("%c", &SHOW_ITERS);
        if( SHOW_ITERS != '\n' )
            while( getchar() != '\n' ) ; /* flush stdin ... */

        if( !(SHOW_ITERS=='1'||SHOW_ITERS=='y'||SHOW_ITERS=='Y') )
        {
            puts("Ok ... 'No' it is ...");
            SHOW_ITERS =0;
        }
        else puts("Ok ... 'Yes' show iterations ... (This may take some time.)");

        printf("\nTo solve a system of linear equations ...\n");
        printf("Enter the number of 'equations/unknowns' to find: ");
        num = 0;
        scanf("%d", &num); /* num of unknowns */
        while( getchar() != '\n' ) ; /* flush stdin ... */

        mat = NULL;
        if( num>0 && (mat = (double**) malloc(num*sizeof(double*))) )
        {
            getInput(num, mat);
            if( checkIfDD(num, mat) ) jacobiCalcDisplay(num,mat);
            else freeMat(num, mat);
        }
        else
        {
            puts("Invalid choice or memory available was exceeded ... Try again.");
            if( mat != NULL ) free( mat );
        }
    }while( more() );
    return 0;
}


void getInput( int numUnKnowns, double** mat )
{
  int i, j;
  FILE *mfile;
  char inputfile[25];
  printf("\nInput file name: ");
  scanf("%s", inputfile);
	mfile = fopen(inputfile, "r");
  if(mfile == NULL){
    printf("Cannot open file %s", inputfile);
    exit(1);
  }
  
	for( i = 0 ; i < numUnKnowns ; i++ )
	{
        mat[i] = (double*) malloc( (numUnKnowns+1)*sizeof(double) );
        puts("");
		for( j = 0 ; j < numUnKnowns+1 ; j++ )
		{
			if( fscanf(mfile, "%lf", &mat[i][j]) != 1 )
            {
                puts("Bad file entry.");
                exit(1);
            }
		}
	}
  
  fclose(mfile);
	printf("\nPress 'Enter' to start iteration ... ");
    getchar();
}

/* Check if the matrix entered is strictly diagonally dominant */
int checkIfDD( int numUnKnowns, double** mat )
{
    int   m, n, dd = 0;
    double* chkdd;
    double* sumdd = (double*) malloc( numUnKnowns*sizeof(double) );
    chkdd = (double*) malloc( numUnKnowns*sizeof(double) );

    for( m = 0 ; m < numUnKnowns ; m++ )
    chkdd[m] = sumdd[m] = 0; /* all set to zero ... */

    printf("\nChecking if the matrix is (strictly) diagonally dominant...");

    for( m = 0 ; m < numUnKnowns ; m++ )
    {
        for( n = 0 ; n < numUnKnowns ; n++ )
            sumdd[m] += fabs(mat[m][n] );
        sumdd[m] -= fabs(mat[m][m]);
        chkdd[m] = fabs(mat[m][m]);

        if(chkdd[m] >= sumdd[m])
        {
            printf("\n%f > %f",chkdd[m],sumdd[m]);
            dd++;
        }
        else printf("\n%f !> %f",chkdd[m],sumdd[m]);
    }

    if(dd == numUnKnowns)
    {
        printf
        (
            "\nYES ..."
            "\nThe matrix is (strictly) diagonally dominant."
            "\nPress 'Enter' to continue with 'Jacobi Iterative Method'...\n"
        );
        getchar();
    }
    else
    {
        printf
        (
            "\nNO ..."
            "\nThe matrix is NOT (strictly) diagonally dominant ... so STOP!"
            "\n(But ... consider exchanging rows in the matrix "
            "and then to try again.)"

            "\n\nPress 'Enter' to continue ... "
        );

        getchar();
        free( sumdd );
        free( chkdd );
        return 0; /* false */
    }

    free( sumdd );
    free( chkdd );
    return 1; /* true */
}

/* uses global SHOW_ITERS ... */
void jacobiCalcDisplay( int numUnKnowns, double** mat )
{
    int* flag;
    int i, j, counter = 0;
    double* res;
    double* var = (double*) malloc( numUnKnowns*sizeof(double) );
    res = (double*) malloc( numUnKnowns*sizeof(double) );
    flag = (int*) malloc( numUnKnowns*sizeof(int) );

    for(i = 0 ; i < numUnKnowns ; i++ )
    var[i] = res[i] = flag[i] = 0;
    printf("The initial value of each array element was set to zero ...\n\n");

    printf( "*********************\n");
    printf( "START CALCULATING ...\n");
    printf( "*********************\n");
    
    clock_t start = clock();

    do
    {
        counter++;
        /* for each iteration keep a copy of the old results ... */
        for(i = 0 ; i < numUnKnowns ; i++ )
        {
            var[i] = res[i];
        }

        if( SHOW_ITERS ) printf("\nIteration number %d ...\n", counter);

        for(i = 0 ; i < numUnKnowns ; i++ ) /* calculation */
        {
            res[i] = mat[i][numUnKnowns];
            for(j = 0 ; j < numUnKnowns ; j++ )
                res[i] = res[i] - mat[i][j]*var[j] ;

            res[i] = res[i] + mat[i][i]*var[i] ;
            res[i] = res[i] / mat[i][i] ;
            if( SHOW_ITERS ) printf("%d = %f\n", i+1, res[i]);
            if( fabs(res[i] - var[i]) < EPSILON ) /* stop condition */
            flag[i]++;

            if( counter==MAX_PRINT_ITERS) SHOW_ITERS = 0;
        }
       
    }while( !checkFlags( numUnKnowns, flag ) );


    printf( "\n********************************\n");
    printf( "The RESULTS of %d ITERATIONS ... \n", counter);
    printf( "********************************\n");

    /*  cross check ...*/

    for( i = 0 ; i < numUnKnowns ; i++)
        var[i] = dotProd( numUnKnowns, mat[i], res );

    showXcheck( numUnKnowns, mat, res, var );

    /* show sol'n vector (again) ... and free up all dynamic memory  */

    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("\nSolution vector ...\n");
    for( i = 0 ; i < numUnKnowns ; i++)
    {
        printf("%d = %+f\n", i+1, res[i]);
        free(mat[i]);
    }
    printf("\nRuntime: %f", time); 
    free( mat );
    free( flag );
    free( res );
    free( var );
}

int checkFlags( int num, int flags[] )
{
    int i;
    for( i=0; i<num; ++ i)
        if( flags[i] == 0 ) return 0;
    return 1;
}

double dotProd( int num, double* v1, double* v2 )
{
    int i;
    double sum =0;
    #pragma acc data copyin(num, v1[0:num-1+1], v2[0:num-1+1]) copy(sum)
    {
    #pragma acc parallel
    {
    #pragma acc loop reduction(+: sum)
    for( i=0; i<num; ++i ) sum += v1[i]*v2[i];
    } // end parallel
    } // end data
    return sum;
}

void showXcheck( int num, double** mat, double* v1, double* v2 )
{
    int i, j;
    puts("\nCross checking ... \nMatrix times sol'n vector ="
         " cal. vector vs. original RHS vector");
    for( i = 0 ; i < num ; i++)
	{
        printf("|");
        for( j =0 ; j < num ; j++ )
            printf("%+9f ", mat[i][j] );
        printf("| |%+9f| |%+9f|vs|%+9f|\n", v1[i], v2[i], mat[i][num]);
    }
}


void freeMat( int num, double** sqMat )
{
    int i;
    for( i=num-1; i>=0; --i )
        free( sqMat[i] );
    free( sqMat );
}