#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>

#define PERCENTILE 1.96

#define PHASE_NUM 101
#define FREQUENCY_NUM 1000


#define START_HZ 0.000125
#define END_HZ 0.1250001
#define INCREMENT_HZ 0.000125

#define RESULT 2

typedef struct sample{
  int gene_count;
  double *sample;
  double sample_mean;
  int timepoint_number;
  double **random_sample;
  double *random_sample_mean;
  int random_sample_number;
  double **zscore;
  double average;
  double *sigma;
}Sample;

int linear_regression(double *x,double *y,unsigned int size,double *a,double *b);
double calculate_SSR(double a,double b,double *x,double *y,int size);
int generate_sample_mean(Sample *samples,int samples_number);
int generate_random_samples(Sample *samples,int samples_number);
int generate_test_cosine(double *test_cosine,double *test_cosine_mean,double *time_points,int sample_num);
int calculate_minimum_SSR(Sample *samples, int sample_num, double *test_cosine, double *test_cosine_mean, int *plotting, double *ans);
void cosopt(double *microarray,double *sigma,double *time_points,int *mrow,int *mcol,int *plotting,double *ans);

int linear_regression(double *x,double *y,unsigned int size,double *a,double *b)
{
  unsigned int i;
  double x_square_total,x_total;
  double y_total,x_multiple_y_total,x_total_square;
  x_square_total = x_total = y_total = x_multiple_y_total = x_total_square = 0;
  for(i=0; i < size;i++){
    x_total += x[i];
    x_square_total += x[i] * x[i];
    y_total += y[i];
    x_multiple_y_total += x[i] * y[i];

  }
  x_total_square = x_total * x_total;


  *a = (size * x_multiple_y_total - x_total * y_total ) / (size * x_square_total - x_total_square);
  *b = (x_square_total * y_total - x_multiple_y_total * x_total ) / (size * x_square_total - x_total_square);
  return 1;
}

/*
 * y =~ ax + b
 */
double calculate_SSR(double a,double b,double *x,double *y,int size){
  int i;
  double distance = 0;
  for(i=0;i<size;i++){
    distance += (y[i] - a * x[i] - b) * (y[i] - a * x[i] - b);
  }
  return distance;
}

int generate_sample_mean(Sample *samples,int samples_number)
{
  int i,j;
  for(i=0;i< samples_number;i++){
    double mean = 0;
    for(j=0;j<samples[i].timepoint_number;j++){
      mean += samples[i].sample[j];
    }
    samples[i].sample_mean = mean / samples[i].timepoint_number;
  }

  return 0;
}

// generate random samples
int generate_random_samples(Sample *samples,int samples_number)
{
  int i,j,k,l;
  double *nr;
  int rt = 0;

  GetRNGstate();

  if(samples != NULL)rt = samples[0].timepoint_number;
  if((rt % 2) == 1){
    rt++;
  }
  nr = (double*)R_alloc(rt,sizeof(double));

  for(i=0;i< samples_number;i++){
    samples[i].random_sample = (double**)R_alloc(samples[i].random_sample_number,sizeof(double*));
    samples[i].random_sample_mean = (double*)R_alloc(samples[i].random_sample_number,sizeof(double));
    for(j=0;j<samples[i].random_sample_number;j++){
      double rm = 0;
      samples[i].random_sample[j] = (double*)R_alloc(samples[i].timepoint_number,sizeof(double));
      memcpy(samples[i].random_sample[j],samples[i].sample,samples[i].timepoint_number*sizeof(double));
      for(k=0;k<samples[i].timepoint_number;k++){
        int n = (int)(unif_rand() * samples[i].timepoint_number);
        double tmp = samples[i].random_sample[j][n];
        samples[i].random_sample[j][n] =  samples[i].random_sample[j][1];
        samples[i].random_sample[j][1] = tmp;
      }
      for(k=0;k<samples[i].timepoint_number;k++){
        //make pseudo gaussian noise
        for(l=0;l<rt;l+=2){
          nr[l] = samples[i].sigma[k] * norm_rand();
          nr[l+1] = samples[i].sigma[k] * norm_rand();
        }
        samples[i].random_sample[j][k] += nr[k];
        rm+=samples[i].random_sample[j][k];
      }
      samples[i].random_sample_mean[j] = (rm / samples[i].timepoint_number);
    }
  }
  PutRNGstate();

  return 0;
}

int generate_test_cosine(double *test_cosine,double *test_cosine_mean,double *time_points,int sample_num)
{
  double i;
  int counter = 0;
  int mean_counter = 0;
  for(i=START_HZ;i<= END_HZ;i+=INCREMENT_HZ){
    int j=0;
    int k=0;

    for(j=50;j >= 1;j--){
      for(k=0;k<sample_num;k++){
        test_cosine[counter] = cos((2*M_PI*i*time_points[k])+(((2.0*M_PI)/100.0) * (double)j));
        test_cosine_mean[mean_counter] += test_cosine[counter];
        counter++;
      }
      test_cosine_mean[mean_counter] = 0;
      mean_counter++;
    }
    //REprintf("%lf:  ",i);
    for(k=0;k<sample_num;k++){
      test_cosine[counter] = cos((2*M_PI*i*time_points[k]));
      //REprintf("time_point=%lf test_cosine=%lf  ,  ",time_points[k],cos((2*M_PI*i*time_points[k])));
      test_cosine_mean[mean_counter] += test_cosine[counter];
      counter++;
    }
    //REprintf("\n");
    test_cosine_mean[mean_counter] = 0;
    mean_counter++;
    for(j=1;j <= 50;j++){
      for(k=0;k<sample_num;k++){
        test_cosine[counter] = cos((2*M_PI*i*time_points[k])-(((2*M_PI)/100) * (double)j));
        test_cosine_mean[mean_counter] += test_cosine[counter];
        counter++;
      }
      test_cosine_mean[mean_counter] = 0;
      mean_counter++;
    }
  }

  return 0;
}

int calculate_minimum_SSR(Sample *samples, int sample_num, double *test_cosine, double *test_cosine_mean, int *plotting, double *ans){
  double i;
  double a,b;
  double distance;
  double *surr_distance=NULL;
  int j,l,m;
  int freq_counter;
  int phase_counter;
  int rsn=0;

  if(samples != NULL){
    rsn = samples[0].random_sample_number;
  }
  if(rsn > 0){
    int x;
    surr_distance = (double*)R_alloc(rsn,sizeof(double));
    for(x=0;x<rsn;x++){
      surr_distance[x] = 0;
    }
  }
  for(l=0;l<sample_num;l++){
    double minimum=100;
    double minimum_freq;
    int minimum_phase;
    freq_counter = 0;
    minimum_freq = minimum_phase = 0;
    for(i=START_HZ;i<= END_HZ;i+=INCREMENT_HZ){
      phase_counter = 0;
      for(j=50;j>=-50;j--){
        double surr_mean=0;
        double surr_sd=0;
        surr_mean = 0;
        surr_sd = 0;
        distance=0;
        linear_regression(samples[l].sample,&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number,&a,&b);
        if( (*plotting==1 && a>=0) || *plotting==0 ){
          distance = calculate_SSR(a,b,samples[l].sample,&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number);
        }else{
          distance = DBL_MAX;
        }
        for(m=0;m<rsn;m++){
          linear_regression(samples[l].random_sample[m],&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number,&a,&b);
          surr_distance[m] = calculate_SSR(a,b,samples[l].random_sample[m],&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number);
          surr_mean += surr_distance[m];
        }
        surr_mean /= rsn;
        for(m=0;m<rsn;m++){
          surr_sd += (surr_distance[m] - surr_mean) * (surr_distance[m] - surr_mean);
        }
        surr_sd = (surr_sd) / (rsn - 1);
        surr_sd = sqrt(surr_sd);
        if((surr_mean - PERCENTILE * surr_sd) > distance && minimum > distance){
          minimum = distance;
          minimum_freq = i;
          minimum_phase = j;
        }
        phase_counter++;
      }
      freq_counter++;
    }
    if(minimum_freq != 0){
      REprintf("minimum_freq=%lf\n",minimum_freq);
      ans[l*RESULT] = minimum_freq;
      ans[l*RESULT+1] = minimum_phase;
    }else{
      REprintf("minimum_freq=0\n");
      ans[l*RESULT] = 0;
      ans[l*RESULT+1] = 0;
    }
  }
  return 0;
}

// arg1: matrix(samples,row=sample,col=time-course)
void cosopt(double *microarray,double *sigma,double *time_points,int *mrow,int *mcol,int *plotting,double *ans){
  double *test_cosine;
  double *test_cosine_mean;
  int i,j;

  int sample_num,gene_num;
  Sample *samples;

  test_cosine = NULL;
  test_cosine_mean = NULL;
  gene_num = *mrow;
  sample_num = *mcol;
  REprintf("gene=%d timecourse=%d\n",gene_num,sample_num);
  REprintf("plotting=%d\n",*plotting);

  test_cosine = (double*)R_alloc(PHASE_NUM * FREQUENCY_NUM * sample_num,sizeof(double));
  test_cosine_mean = (double*)R_alloc(PHASE_NUM * FREQUENCY_NUM,sizeof(double));
  //warning("number of test cosine = %d\n",PHASE_NUM * FREQUENCY_NUM * sample_num);
  if(test_cosine == NULL || test_cosine_mean == NULL){
    error("Memory allocation error\n");
  }
  for(i=0;i<PHASE_NUM * FREQUENCY_NUM;i++){
    test_cosine_mean[i] = 0;
  }
  generate_test_cosine(test_cosine,test_cosine_mean,time_points,sample_num);

  samples = (Sample*)R_alloc(gene_num,sizeof(Sample));
  for(i=0;i<gene_num;i++){
    samples[i].gene_count = i;
    samples[i].sample = (double*)R_alloc(sample_num,sizeof(double));
    samples[i].sigma = (double*)R_alloc(sample_num,sizeof(double));
    for(j=0;j<sample_num;j++){
      samples[i].sample[j] = microarray[i+(gene_num)*j];
      samples[i].sigma[j] = sigma[i+(gene_num)*j];
    }
    samples[i].average = 0;
    samples[i].random_sample_number = 1000;
    samples[i].timepoint_number = sample_num;
    samples[i].random_sample = NULL;
  }
  generate_sample_mean(samples,gene_num);
  generate_random_samples(samples,gene_num);
  calculate_minimum_SSR(samples,gene_num,test_cosine,test_cosine_mean,plotting,ans);
}
