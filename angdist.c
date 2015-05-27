/* angdist.c - calcolo della distribuzione della distanza angolare tra coppie di punti nello spazio.
 *
 *
 * There are no restrictions about the number of threads that can be used
 *
 * Compile with:
 * mpicc -Wall -lm -g angdist.c -o angdist
 *
 * Run with:
 * mpirun -n <number_of_threads> angdist <number_of_intervals> <number_of_stars> <input_file>
 *
 * Author: Ettore Di Giacinto <mudler@sabayon.org>
 * 
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#define mpi_root 0
#define GNUPLOT "gnuplot -persist"

typedef struct star
{
    double x;
    double y;
    double z;
} star_t;

void *load_tree(star_t *, unsigned long, char[]);
int file_exist (char *);
double star_angular_distance(star_t , star_t );
void draw_histogram( unsigned int *, int size);


/*===============================================================*/
/* load_tree - load the file and fill up the memory              */
/* @args:   stars_t  stars pointer, char[] the file name     */
/* returns: stars_t pointer with all data filled               */
/*===============================================================*/

void *load_tree(star_t *stars, unsigned long stars_n, char file[])
{
    int scanned_stars = 0;
    FILE *file_in;  // pointer to filehandle
    file_in = fopen(file, "r");   //open file
    stars = malloc(sizeof(star_t) * stars_n); //And prepare the memory for the stars
    int x;
        for (x=0; x< stars_n; x++){
             if ( fscanf(
                    file_in,
                    "%lf %lf %lf",
                    &stars[scanned_stars].x ,
                    &stars[scanned_stars].y,
                    &stars[scanned_stars].z) != EOF)
        {
            scanned_stars++;// If the fscanf goes well we can keep going
        }
        }

    if(scanned_stars != stars_n)
    {
        printf("[!] !! Something really wrong happened, your dataset could be malformed, or maybe it's just me !!");
    }
    fclose(file_in);
    return (stars);
}




/*===============================================================*/
/* file_exist - if the reading of attribute it's ok              */
/*              the file is going to exist                       */
/* @args:   char[] file name                                     */
/* returns: int 1 if it exists, 0 otherwise                      */
/*===============================================================*/

int file_exist (char *filename)
{
    struct stat   buffer;
    return (stat (filename, &buffer) == 0);  //The stat returns 0 if the attributes could be read, -1 otherwise
}

/*===============================================================*/
/* star_angular_distance - computes the angular distance between */
/*                         two stars                             */
/* @args:   star_t star_i, star_t star_j                         */
/* returns: angular distance (double)                            */
/*===============================================================*/

double star_angular_distance(star_t star_i, star_t star_j)
{
    return acos(
               ((star_i.x * star_j.x) + (star_i.y * star_j.y ) + (star_i.z * star_j.z ))
               /
               sqrt(
                // pow() :  gives more precision, but it's too resource-expensive
                   //        (pow(  star_i.x, 2) + pow( star_i.y , 2)
                   //        + pow( star_i.z, 2)  ) *
                   //      (pow(  star_j.x, 2) + pow( star_j.y , 2)
                   //      + pow( star_j.z, 2)  )     )
                   (  star_i.x *  star_i.x + star_i.y * star_i.y +  star_i.z *  star_i.z)
                   *
                   (star_j.x * star_j.x + star_j.y * star_j.y + star_j.z * star_j.z)
               )
           );
}

/*===============================================================*/
/* draw_histogram - opens a pipe to gnuplot, print the histogram */
/*                  and save all the data to                     */
/*		    star_angular_distance_distribution.dat       */
/* @args:   unsigned int *histogram, int histogram_size          */
/* returns: nothing                                              */
/*===============================================================*/

void draw_histogram( unsigned int *histogram, int size)
{
    int status;
    FILE *gp;

    gp = popen(GNUPLOT, "w"); /* 'gp' is the pipe descriptor */
    if (gp == NULL)
    {
        printf("[!] Error opening pipe to GNU plot. Check if you have it! \n");
    }

    FILE *f = fopen("star_angular_distance_distribution.dat", "w");

    if (f == NULL)
    {
        printf("Error opening file star_angular_distance_distribution.dat!\n");
    }

    int i;
    for(i = 0; i < size; i++)
    {
        fprintf(f, "%f %u\n", (double) ((double) i / size * M_PI), histogram[i]);
    }
    fclose(f);

    /* GNUPLOT OPTIONS */

    fprintf(gp,"set term png truecolor\n");
    fprintf(gp, "set output 'star_angular_distance_distribution.png'\n");
    fprintf(gp,"set xlabel 'Distance'\n");
    fprintf(gp,"set ylabel 'Angular distances'\n");
    fprintf(gp, "set xrange [0:%f]\n", M_PI);
    fprintf(gp,"set grid\n");
    fprintf(gp,"set boxwidth 0.95 relative\n");
    fprintf(gp,"set style fill transparent solid 0.5 noborder\n");
    fprintf(gp, "plot 'star_angular_distance_distribution.dat' with boxes lc rgb 'blue' notitle\n");


    fprintf(gp,"exit\n");
    status = pclose(gp);

    if (status == -1)
    {
        printf("[!] Something went wrong using gnuplot, i'm sorry.\n");
    }
}


int main(int argc, char *argv[])
{

    int intervals;
    unsigned long stars_n;
    unsigned int *histogram = NULL, *histogram_local = NULL;
    char *file_name;

    if ( 4 == argc )
    {
        intervals = atoi(argv[1]);
        stars_n = atoi(argv[2]);
        file_name = argv[3];
    }
    if ( !file_exist( file_name ) ) //let's check the existence of the file first
    {

        printf("/!\\ Can't access to %s\n", file_name);
        return 2;

    }

    /*===============================================================*/
    /* MPI initialization                                            */
    /*===============================================================*/

    MPI_Init(&argc, &argv);

    /*===============================================================*/
    /* Vars for MPI threads                                          */
    /*===============================================================*/

    unsigned long             i, j;

    int taskid, numtasks,
        proclen;
    char            procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(procname, &proclen);

    int             blockcounts[1];
    MPI_Aint        offsets[1];
    MPI_Datatype     startype, oldtypes[1];
    unsigned long total_distances = (stars_n * (stars_n - 1)) / 2; // N*(N-1)/2
    unsigned long             distances_to_compute = (total_distances / numtasks) + (taskid < total_distances % numtasks); // distributing equally the remaining tasks
    double             interval_size = M_PI / (double) intervals ;
    double          starttime, endtime, total_time;
    star_t     *stars = NULL;
    unsigned long localstart = 0, localstop = 0; 

    /*===============================================================*/
    /* Definition and commit of MPI structures                       */
    /*===============================================================*/


    /* computation of the 3 MPI_DOUBLE fields x, y, z byte occupation for stars */
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;  // on first block(0) there are 3 MPI_DOUBLE doubles
    blockcounts[0] = 3;
    MPI_Type_struct(1, blockcounts, offsets, oldtypes, &startype); /* 1 is the count of the types on arraytypes */
    MPI_Type_commit(&startype);

    total_time = MPI_Wtime();


    /*===============================================================*/


    stars = (star_t *) malloc( stars_n * sizeof(star_t));


    if( taskid == mpi_root)
    {
        /* Load of stars on "stars" */
        stars = load_tree(stars, stars_n, file_name);
        printf("[**] Intervals: %d\n[**] Number of Stars: %lu\n[**] Dataset: %s\n", intervals, stars_n, file_name);
        printf("[**] Total distances: %lu\n", total_distances);
        histogram = calloc(intervals, sizeof(unsigned int));
    }
    else
    {
        localstart = (taskid ) * ((total_distances / (numtasks)) + (taskid - 1 < total_distances % (numtasks)) * 2) ;
    }
    histogram_local = calloc(intervals, sizeof(unsigned int));

    localstop = localstart + distances_to_compute;


    printf("-[ Process:  %d ]- Distance start: %lu \t Distance ending %lu -----  %lu to compute\n", taskid, localstart, localstop, distances_to_compute);

    /*===============================================================*/
    /* Broadcasting all the stars_n stars on every process           */
    /*===============================================================*/

    MPI_Bcast(stars,  stars_n, startype, mpi_root, MPI_COMM_WORLD);

    starttime = MPI_Wtime(); 

    long double M = (stars_n - 1); //avoid diagonal

    for (i = localstart; i < localstop; ++i)
    {
        long double gt_i = M * (M + 1) / 2 - 1 - i;
        long double K = floor((sqrt(8 * gt_i + 1) - 1) / 2);

        double distance = star_angular_distance(

                              stars[(unsigned long) (M - 1 - K)],// Resolves the star in the row
                              stars[
                                  (unsigned long)
                                  (
                                      i - M * (M + 1) / 2 + (K + 1) * (K + 2) / 2
                                  )
                              ] // Resolves the star in the column

                          );
        unsigned int index = floor( distance / interval_size); // find the corresponding interval
        histogram_local[index]++;
    }

    endtime   = MPI_Wtime();
    MPI_Barrier( MPI_COMM_WORLD ); // So output will be syncronized

    printf("-[ Process:  #%d ]- Operation took: %f.3 s \n\n", taskid, endtime - starttime);

    starttime = MPI_Wtime();

    MPI_Reduce(histogram_local,
               histogram,
               intervals,
               MPI_UNSIGNED,
               MPI_SUM,
               mpi_root,
               MPI_COMM_WORLD);

    endtime   = MPI_Wtime();
    printf("-[ Process: #%d ]- finished calculating in %f seconds\n", taskid, endtime - starttime);


    if (taskid == mpi_root)
    {
        unsigned long sum = 0;
        for (j = 0; j < intervals; ++j)
        {
            sum += histogram[j];

            printf("F_%lu\t=\t%u\n", j, histogram[j]);
        }
        draw_histogram(histogram, intervals);
        free(histogram);
        total_time = MPI_Wtime() - total_time;
        printf("-[ Stats ]- Total Execution time is of %fs\n", total_time);
        if(sum != total_distances)
        {
            // Just a minimal checking if we have done something really wrong
            printf("...lies... LIES...... LIES!!\n");
        }
    }

    free(histogram_local);
    free(stars);

    MPI_Finalize();
    return 0;
}
