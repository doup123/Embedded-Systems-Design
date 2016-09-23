/*Parallel Hierarchical One-Dimensional Search for motion estimation*/
/*Initial algorithm - Used for simulation and profiling             */

#include <time.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#define MILLION 1000000L
#define N 144     /*Frame dimension for QCIF format*/
#define M 176     /*Frame dimension for QCIF format*/
#define p 7       /*Search space. Restricted in a [-p,p] region around the
                    original location of the block.*/
int B_x,B_y;//Block dimensions
void read_sequence(int current[144][176], int previous[144][176])
{ 
	FILE *picture0, *picture1;
	int i, j;

	if((picture0=fopen("akiyo0.y","rb")) == NULL)
	{
		printf("Previous frame doesn't exist.\n");
		exit(-1);
	}

	if((picture1=fopen("akiyo1.y","rb")) == NULL) 
	{
		printf("Current frame doesn't exist.\n");
		exit(-1);
	}

  /*Input for the previous frame*/
  for(i=0; i<144; i++)
  {
    for(j=0; j<176; j++)
    {
      previous[i][j] = fgetc(picture0);
      current[i][j] = fgetc(picture1);
    }
  }



	fclose(picture0);
	fclose(picture1);
}


void phods_motion_estimation(int current[144][176], int previous[144][176],
    int **vectors_x,int **vectors_y)
{
int  p1, p2, q2, distx, disty,disty0,distx0,distx_,disty_, min1,min2,bestx, besty;
register int x,y,k,l,v1,v2;
 distx = 0;
 disty = 0;
 
  for(x=0; x<N/B_x; x++)
  {
    for(y=0; y<M/B_y; y++)
    {
		register int Bx=B_x*x;
		register int By=B_y*y;
		v1 = vectors_x[x][y] +Bx;
		v2 = vectors_y[x][y] +By;
        min1 = 255*B_x*B_y;
        min2 = 255*B_x*B_y;
       //////////////////////////////////
	//S = 4;//S=-4,0,4
       //////////////////////////////////

         distx=0;
		disty=0;
		disty0=0;
		distx0=0;
		disty_=0;
		distx_=0;
        for(k=0; k<B_x; k++)     
          {
            for(l=0; l<B_y; l++)
            {
              p1 = current[Bx+k][By+l];

              if((v1 -4 + k) < 0 ||
                  (v1 -4 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1-4+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 -4 + l) < 0 ||
                  (v2 -4 + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
	      {
                q2 = previous[v1+k][v2-4+l];
              }
              disty_ += abs(p1-q2);
              distx_ += abs(p1-p2);
              if((v1  + k) < 0 ||
                  (v1 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2  + l) < 0 ||
                  (v2  + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
	      {
                q2 = previous[v1+k][v2+l];
              }
	       disty0 += abs(p1-q2);
               distx0 += abs(p1-p2);
               if((v1 +4 + k) < 0 ||
                  (v1 + 4 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1+4+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 +4 + l) < 0 ||
                  (v2 +4 + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
 	      {
                q2 = previous[v1+k][v2+4+l];
              }
	      disty += abs(p1-q2);
              distx += abs(p1-p2);
	   }
          }
	/*find the minimum of disty,disty0,disty_ and the same in x */	  
	if (disty<disty0)
	{
		if (disty<disty_)
		{
			besty = 4;		
		}		
		else
		{
			besty = -4;
		}
	}
	else
	{
		if(disty0 < disty_)
		{
			besty = 0;
		}
		else
		{
			besty = -4;
		}
	}
	if (distx<distx0)
	{
		if (distx<distx_)
		{
			bestx = 4;		
		}		
		else
		{
			bestx = -4;
		}
	}
	else
	{
		if(distx0 < distx_)
		{
			bestx = 0;
		}
		else
		{
			bestx = -4;
		}
	}
		/*Add each minimum amount to the vectors*/
		v1 += bestx;
        v2 += besty;
       
       //////////////////////////////////
      //  S = 2;//S=-2,0,2
       //////////////////////////////////
		bestx=0;
		besty=0;	
		disty0=0;
		distx0=0;
		distx=0;
		disty=0;
		disty_=0;
		distx_=0;
        for(k=0; k<B_x; k++)     
          {
            for(l=0; l<B_y; l++)
            {
              p1 = current[Bx+k][By+l];

              if((v1  -2 + k) < 0 ||
                  (v1-2 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1-2+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 -2 + l) < 0 ||
                  (v2 -2 + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
  	      {
                q2 = previous[v1+k][v2-2+l];
              }

              disty_ += abs(p1-q2);
              distx_ += abs(p1-p2);
	         if((v1 + 2 + k) < 0 ||
                  (v1 + 2 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[Bx+v1+2+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 +2 + l) < 0 ||
                  (v2 +2 + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
	      {
                q2 = previous[v1+k][v2+2+l];
              }

              disty += abs(p1-q2);
              distx += abs(p1-p2);
			if((v1  + k) < 0 ||
                  (v1 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2  + l) < 0 ||
                  (v2  + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
	      {
                q2 = previous[v1+k][v2+l];
              }
	       disty0 += abs(p1-q2);
               distx0 += abs(p1-p2);
	      

            }
          }
	/*find the minimum of disty,disty0,disty_ and the same in x */	  
        if (disty<disty0)
	{
		if (disty<disty_)
		{
			besty = 2;		
		}		
		else
		{
			besty = -2;
		}
	}
	else
	{
		if(disty0 < disty_)
		{
			besty = 0;
		}
		else
		{
			besty = -2;
		}
	}
	if (distx<distx0)
	{
		if (distx<distx_)
		{
			bestx = 2;		
		}		
		else
		{
			bestx = -2;
		}
	}
	else
	{
		if(distx0 < distx_)
		{
			bestx = 0;
		}
		else
		{
			bestx = -2;
		}
	}
		/*Add each minimum amount to the vectors*/
        v1 += bestx;
        v2 += besty;

       //////////////////////////////////
      //  S = 1;//S=-1,0,1
       //////////////////////////////////
     bestx=0;
	besty=0;
	disty0=0;
	distx0=0;	
	min1=distx0;
        min2=disty0;
        distx=0;
	disty=0;
	disty_=0;
	distx_=0;
          
          for(k=0; k<B_x; k++)     
          {
            for(l=0; l<B_y; l++)
            {
              p1 = current[Bx+k][By+l];

              if((v1 -1 + k) < 0 ||
                  (v1 -1 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1-1+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 -1 + l) < 0 ||
                  (v2 -1 + l) > (176-1))
              {
                q2 = 0;
              } 	
	      else 
	      {
                q2 = previous[v1+k][v2-1+l];
              }

              disty_ += abs(p1-q2);
              distx_ += abs(p1-p2);
              if((v1 + 1 + k) < 0 ||
                  (v1 + 1 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1+1+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2 +1 + l) < 0 ||
                  (v2 +1 + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
	      {
                q2 = previous[v1+k][v2+1+l];
              }

              disty += abs(p1-q2);
              distx += abs(p1-p2);
		 if((v1  + k) < 0 ||
                  (v1 + k) > (144-1) ||
                  (v2 + l) < 0 ||
                  (v2 + l) > (176-1))
              {
                p2 = 0;
              } 
	      else 
	      {
                p2 = previous[v1+k][v2+l];
              }
              if((v1 + k) <0 ||
                  (v1 + k) > (144-1) ||
                  (v2  + l) < 0 ||
                  (v2  + l) > (176-1))
              {
                q2 = 0;
              } 
	      else 
 	      {
                q2 = previous[v1+k][v2+l];
              }
	       disty0 += abs(p1-q2);
               distx0 += abs(p1-p2);
            }
          }
 	/*find the minimum of disty,disty0,disty_ and the same in x */	  
	if (disty<disty0)
	{
		if (disty<disty_)
		{
			besty = 1;		
		}		
		else
		{
			besty = -1;
		}
	}
	else
	{
		if(disty0 < disty_)
		{
			besty = 0;
		}
		else
		{
			besty = -1;
		}
	}
	if (distx<distx0)
	{
		if (distx<distx_)
		{
			bestx = 1;		
		}		
		else
		{
			bestx = -1;
		}
	}
	else
	{
		if(distx0 < distx_)
		{
			bestx = 0;
		}
		else
		{
			bestx = -1;
		}
	}
	/*Add each minimum amount to the vectors*/
	v1 += bestx;
    v2 += besty; 
/*keep only the value of the vectors*/
	vectors_x[x][y] = v1 - Bx;
	vectors_y[x][y] = v2 -By;
    }
  }
} 

void usage(char *argv0)
{
	fprintf(stderr, "Usage: %s \n\n"
		"Exactly two arguments required:\n"
		"Bx: the horizontal dimension of the block size, one of the dividers of 144\n"
		"By: the vertical dimension of the block size, one of the dividers of 176\n",
		argv0);
	exit(1);
}

int main(int argc, char *argv[])
{
//int B;  

int current[144][176], previous[144][176], i, j;
	if (argc != 3){
		/*Incorrect number of arguments*/
		usage(argv[0]);
	}
	/*Get the arguments*/
	B_x = atoi(argv[1]);
	B_y = atoi(argv[2]);
	/*The arguments need to be proper dividers of 144 or 176*/
	if ( (144 % B_x) != 0) {
		fprintf(stderr, "`%s' is not valid for horizontal dimension size\n", argv[1]);
		exit(1);
	}
	if((176 % B_y) != 0){
		fprintf(stderr, "`%s' is not valid for vertical dimension size\n", argv[2]);
		exit(1);
	}
	int *arr1[144/B_x],*arr2[144/B_x] ;
	/*Allocate memory space*/
	for(i=0;i<(144/B_x);i++)
	{
		arr1[i] = (int*)malloc((176/B_y)*sizeof(int));
		arr2[i] = (int*)malloc((176/B_y)*sizeof(int));
	}
	/*Initialize vectors*/
	for(i = 0; i< (144/B_x); i++){
		for(j=0; j<(176/B_y); j++){
			arr1[i][j] = 0;
		}
	}  
    struct timeval tts,ttf;
	double ttotal;
	read_sequence(current,previous);
	gettimeofday(&tts,NULL);
  	phods_motion_estimation(current,previous,arr1,arr2);
	gettimeofday(&ttf,NULL);
	/*total time needed*/
	ttotal= MILLION* (ttf.tv_sec - tts.tv_sec)+ ttf.tv_usec-tts.tv_usec;
	int ttotal1;
	ttotal1 = (int)ttotal;
	printf("%d",ttotal1);
	/*free the allocated memory*/
for ( i = 0; i < 144/B_x; i++)
{
    int* currentIntPtr = arr1[i];
    free(currentIntPtr);
    currentIntPtr = arr2[i];
	free(currentIntPtr);
}
  return 0;
}
