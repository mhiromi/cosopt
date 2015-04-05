#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <dSFMT.h>
#include <SFMT.h>
#include <math.h>

#define SIGMA 0.3

#define MAX 1024
#define MAXGENE 40960

#define PERCENTILE 1.96

#define PHASE_NUM 101
#define FREQUENCY_NUM 1000

#define START_HZ 0.000125
#define END_HZ 0.1250001
#define INCREMENT_HZ 0.000125

double global_sigma = SIGMA;

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
  double sigma;
}Sample;

int linear_regression(double *x,double *y,unsigned int size,double *a,double *b)
{
  int i,j,k;
  double x_square_total,x_total;
  double y_total,x_multiple_y_total,x_total_square;
  x_square_total = x_total = y_total = x_multiple_y_total = x_total_square = 0;
  for(i=0;i<size;i++){
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

int check_expression_file(char *filename,int *gene_num,int *sample_num)
{
  FILE *fp;
  char line[MAX];
  int i;
  int gene_count = 0;
  int first_line_flag =0;
  char microarray_name[MAX];
  if((fp = fopen(filename,"r")) < 0){
    fprintf(stderr,"Can't open %s",filename);
    exit(1);
  }
  while (fgets(line, MAX, fp) != NULL) {
    int cnum = strlen(line);
    int word_count = 0;
    char word[MAX];
    int sn = 0;
    int cc = 0;
    memset(word,0,MAX);
    if(first_line_flag == 0){
      first_line_flag++;
    }else{
      for(i=0;i < cnum;i++){
        if(line[i] == ',' || line[i] == '\n' || line[i] == '\t' || line[i] == '#'){
          word[cc] = '\0';
          if(word_count == 0){
            strcpy(microarray_name,word);
          }else{
            sn++;
          }
          word_count++;
					if(line[i] == '#') break;
				}else if(line[i] == ' '){
					continue;
        }else{
          char c = line[i];
          word[cc] = c;
        }
      }
      if(sn > *sample_num)*sample_num = sn;
      gene_count++;
    }
  }
  *gene_num = gene_count;
  return(0);
}

int load_expression_file(char *filename, double* microarray, char** microarray_name,double *time_points,int gene_num,int sample_num)
{
  FILE *fp;
  char line[MAX];
  int i;
  int gene_count = 0;
  int firstline_flag = 0;
	int comment_flag = 0;
  char word[MAX];
  if((fp = fopen(filename,"r")) < 0){
    fprintf(stderr,"Can't open %s",filename);
    exit(1);
  }
  while (fgets(line, MAX, fp) != NULL) {
    if(firstline_flag == 0){
      int cnum = strlen(line);
      int t=0;
      int cc=0;
      memset(word,0,MAX);
      for(i=0;i < cnum;i++){
        if(line[i] == ',' || line[i] == '\n' || line[i] == '\t'){
          word[cc] = '\0';
          time_points[t] = atof(word);
          t++;
          memset(word,0,MAX);
          cc=0;
        }else{
          char c = line[i];
          word[cc] = c;
          cc++;
        }
      }
      firstline_flag++;
    }else{
      int cnum = strlen(line);
      int word_count = 0;
      int sn = 0;
      int cc=0;
      memset(word,0,MAX);
      for(i=0;i < cnum;i++){
        if(line[i] == ',' || line[i] == '\n' || line[i] == '\t' || line[i] == '#'){
          word[cc] = '\0';
          if(word_count == 0){
            strcpy(microarray_name[gene_count],word);
          }else{
            double tmp;
            tmp = atof(word);
            microarray[(gene_count*sample_num)+sn] = atof(word);
            sn++;
          }
          memset(word,0,MAX);
          word_count++;
          cc=0;
					if(line[i] == '#') break;
				}else if(line[i] == ' '){
					continue;
				}else{
          char c = line[i];
          word[cc] = c;
          cc++;
        }
      }
      gene_count++;
    }
  }
  return(0);
}

int generate_sample_mean(Sample *samples,int samples_number)
{
  int i,j,k,l;
  for(i=0;i< samples_number;i++){
    double mean = 0;
    for(j=0;j<samples[i].timepoint_number;j++){
      mean += samples[i].sample[j];
    }
    samples[i].sample_mean = mean / samples[i].timepoint_number;
    //printf("%f\n",samples[i].sample_mean);
  }
}
// CORRCOSのやり方でランダムのサンプルを生成する
// タイムポイントはどのサンプルも同じだと仮定する
int generate_random_samples(Sample *samples,int samples_number)
{
  int i,j,k,l;
  double *nr;
  int rt = 0;
  int_init_gen_rand(time(NULL));
  init_gen_rand(time(NULL));

  if(samples != NULL)rt = samples[0].timepoint_number;
  if((rt % 2) == 1){
    rt++;
  }
  nr = (double*)malloc(rt * sizeof(double));

  for(i=0;i< samples_number;i++){
    samples[i].random_sample = (double**)malloc(samples[i].random_sample_number * sizeof(double*));
    samples[i].random_sample_mean = (double*)malloc(samples[i].random_sample_number * sizeof(double));
    for(j=0;j<samples[i].random_sample_number;j++){
      double rm = 0;
      samples[i].random_sample[j] = (double*)malloc(samples[i].timepoint_number * sizeof(double));
      memcpy(samples[i].random_sample[j],samples[i].sample,samples[i].timepoint_number*sizeof(double));
      // タイムポイントの数だけとにかくシャッフリング
      for(k=0;k<samples[i].timepoint_number;k++){
        int n = gen_rand32() % samples[i].timepoint_number;
        double tmp = samples[i].random_sample[j][n];
        samples[i].random_sample[j][n] =  samples[i].random_sample[j][1];
        samples[i].random_sample[j][1] = tmp;
      }
      // pseudo-gaussianノイズを入れる
      for(k=0;k<samples[i].timepoint_number;k++){
        double r1,r2;
        //Box-Muller transform
        for(l=0;l<rt;l+=2){
          r1 = genrand_open_close();
          r2 = genrand_open_close();
          // できた乱数をaverageとsigmaで変換しておく
          nr[l] = samples[i].sigma * sqrt(-2*log(r1))*sin(2*M_PI*r2) + samples[i].average;
          nr[l+1] = samples[i].sigma * sqrt(-2*log(r1))*cos(2*M_PI*r2) + samples[i].average;
        }
        samples[i].random_sample[j][k] += nr[k];
        rm+=samples[i].random_sample[j][k];
        //printf("%f\n", samples[i].random_sample[j][k]);
      }
      samples[i].random_sample_mean[j] = (rm / samples[i].timepoint_number);
    }
  }
  free(nr);
}

int generate_test_cosine(double *test_cosine,double *test_cosine_mean,double *time_points,int sample_num)
{
  double i;
  int counter = 0;
  int mean_counter = 0;
  for(i=START_HZ;i<= END_HZ;i+=INCREMENT_HZ){
    int j=0;
    int k=0;

    //printf("%f\n",i);
    for(j=50;j >= 1;j--){
      for(k=0;k<sample_num;k++){
        test_cosine[counter] = cos((2*M_PI*i*time_points[k])+(((2*M_PI)/100) * (double)j));
        test_cosine_mean[mean_counter] += test_cosine[counter];
        counter++;
        //printf("%lf %d\n",time_points[k],k);
      }
      test_cosine_mean[mean_counter] = 0;
      mean_counter++;
    }
    for(k=0;k<sample_num;k++){
      test_cosine[counter] = cos((2*M_PI*i*time_points[k]));
      test_cosine_mean[mean_counter] += test_cosine[counter];
      counter++;
    }
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
}

int calculate_minimum_SSR(Sample *samples, int sample_num, double *test_cosine, double *test_cosine_mean, char **microarray_name){
  double i;
  double a,b;
  double distance;
  double *surr_distance;
  int j,k,l,m;
  int counter;
  int freq_counter;
  int phase_counter;
  int rsn;
  if(samples != NULL)rsn = samples[0].random_sample_number;
  if(rsn > 0){
    int x;
    surr_distance = (double*)malloc(rsn * sizeof(double));
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
        linear_regression(samples[l].sample,&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number,&a,&b);
        distance = calculate_SSR(a,b,samples[l].sample,&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number);
        //				fprintf(stderr,"[%lf][%d] %lf %lf %lf\n",1.0/i,j,a,b,distance);
        for(m=0;m<rsn;m++){
          linear_regression(samples[l].random_sample[m],&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number,&a,&b);
          surr_distance[m] = calculate_SSR(a,b,samples[l].random_sample[m],&test_cosine[(freq_counter * 101 * samples[l].timepoint_number) + (phase_counter * samples[l].timepoint_number)],samples[l].timepoint_number);
          surr_mean += surr_distance[m];
          //				  if(m == 0)fprintf(stderr,"[random][%lf][%d]%lf %lf %lf\n",1.0/i,j,a,b,surr_distance[m]);
        }
        surr_mean /= rsn;
        for(m=0;m<rsn;m++){
          surr_sd += (surr_distance[m] - surr_mean) * (surr_distance[m] - surr_mean);
        }
        surr_sd = (surr_sd) / (rsn - 1);
        surr_sd = sqrt(surr_sd);
        //				fprintf(stderr,"%lf\n",surr_sd);
        //				fprintf(stderr,"[0.95][%lf][%d]%lf %lf %lf %lf\n",1.0/i,j,distance,surr_mean,PERCENTILE * surr_sd,(surr_mean - PERCENTILE * surr_sd));
        if((surr_mean - PERCENTILE * surr_sd) > distance && minimum > distance){
          minimum = distance;
          minimum_freq = i;
          minimum_phase = j;
          //					fprintf(stderr,"[mini][%lf][%d] %lf\n",1.0/(minimum_freq),minimum_phase,minimum);
        }
        phase_counter++;
      }
      freq_counter++;
    }
    if(minimum_freq != 0){
      //		  fprintf(stderr,"%s,%lf,%d,%lf\n",microarray_name[samples[l].gene_count],1.0/(minimum_freq),minimum_phase,minimum);
      fprintf(stdout,"%s,%lf,%d,%lf\n",microarray_name[samples[l].gene_count],1.0/(minimum_freq),minimum_phase,minimum);
    }
  }
}

int calculate_zscore(Sample *samples,int sample_num,double *test_cosine,double *test_cosine_mean,char **microarray_name){
  double i;
  int j,k,l,m,n;
  int counter;
  int mean_counter;
  double *csurr,*csurr_autocorrelation;
  int rsn=0;
  FILE *fp;

  char filename[MAX];

  if(samples != NULL)rsn = samples[0].random_sample_number;
  if(rsn > 0){
    int x;
    csurr = (double*)malloc(rsn * sizeof(double));
    csurr_autocorrelation = (double*)malloc(rsn * sizeof(double));
    for(x=0;x<rsn;x++){
      csurr[x] = csurr_autocorrelation[x] = 0;
    }
  }
  //printf("random_sample_number = %d\n",rsn);
  for(l=0;l<sample_num;l++){
    sprintf(filename,"%s.txt",microarray_name[samples[l].gene_count]);
    if((fp = fopen(filename,"w")) < 0){
      fprintf(stderr,"Can't open %s",filename);
      exit(1);
    }

    counter = 0;
    mean_counter = 0;
    int freq_counter = 0;
    for(i=START_HZ;i<= END_HZ;i+=INCREMENT_HZ){
      int j=0;
      int k=0;
      int x=0;
      int phase_counter=0;
      //		printf("%f\n",i);

      for(j=50;j >= 1;j--){
        double csurr_total=0;
        double csurr_total_sum=0;
        double ccorig = 0;
        double sample_autocorrelation=0;
        double cosine_autocorrelation=0;
        for(x=0;x<rsn;x++){
          csurr_autocorrelation[x] = csurr[x] = 0;
        }
        for(k=0;k<samples[l].timepoint_number;k++){
          ccorig += (samples[l].sample[k] - samples[l].sample_mean) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          sample_autocorrelation += (samples[l].sample[k] - samples[l].sample_mean) * (samples[l].sample[k] - samples[l].sample_mean);
          cosine_autocorrelation += (test_cosine[counter] - test_cosine_mean[mean_counter]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]test_cosine=%f cosine_mean=%f sample=%f sample_mean=%f\n",i,j,k,test_cosine[counter] ,test_cosine_mean[mean_counter],samples[l].sample[k],samples[l].sample_mean);
          //fprintf(stderr,"[%lf][%d][%d]test_cosine=%f cosine_mean=%f\n",i,j,k,test_cosine[counter] ,test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]     sample=%f sample_mean=%f\n",i,j,k,samples[l].sample[k] , samples[l].sample_mean);
          //printf("value[%d][%d]=%lf\n",l,k,samples[l].sample[k]);
          for(m=0;m<rsn;m++){
            csurr[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
            //fprintf(stderr,"csurr[%d]=%lf\n",m,(samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]));
            csurr_autocorrelation[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]);
            //printf("[%d] sample=%lf sample_mean=%lf test_cosine=%lf cosine_mean=%lf csurr=%lf\n",m,samples[l].random_sample[m][k] ,samples[l].sample_mean,test_cosine[counter],test_cosine_mean[counter],csurr[m]);
          }
          //test_cosine[counter] = cos((M_PI*i*time_points[k])+(((M_PI*i)/100) * (double)j));
          counter++;
        }
        ccorig /= ((double)samples[l].timepoint_number - 1);
        ccorig /= (sqrt((sample_autocorrelation) / ((double)samples[l].timepoint_number - 1)) * sqrt((cosine_autocorrelation) / ((double)samples[l].timepoint_number - 1)));
        //			fprintf(stderr,"ccorig=%f\n",ccorig);
        for(m=0;m<rsn;m++){
          csurr[m] /= ((double)samples[l].timepoint_number - 1);
          csurr[m] /= (sqrt(csurr_autocorrelation[m] / ((double)samples[l].timepoint_number - 1)) * sqrt(cosine_autocorrelation / ((double)samples[l].timepoint_number - 1)));
          //fprintf(stderr,"[%lf][%d]csurr_auto=%lf cosine_auto=%lf ccorig=%lf csurr[%d]=%lf\n",i,j,sqrt(csurr_autocorrelation[m]),sqrt(cosine_autocorrelation),m,ccorig,csurr[m]);
          csurr_total += csurr[m];
        }
        csurr_total /= (double)samples[l].random_sample_number;
        for(m=0;m<rsn;m++){
          csurr_total_sum += (csurr[m] - csurr_total) * (csurr[m] - csurr_total);
          //	fprintf(stderr,"sd[%d]=%f csurr_total_sum=%f\n",m,(csurr[m] - csurr_total) * (csurr[m] - csurr_total),csurr_total_sum);
        }
        csurr_total_sum /= (double)(samples[l].random_sample_number-1);
        csurr_total_sum = sqrt(csurr_total_sum);
        //			fprintf(stderr,"ccorig=%lf csurr_total[%d]=%lf SDsurr=%lf\n",ccorig,j,csurr_total,csurr_total_sum);
        //			fprintf(stderr,"[%d][%lf][%d] test_zscore = %lf\n",samples[l].gene_count,1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);

        //			fprintf(fp,"[%d][%lf][%d] test_zscore = %lf\n",samples[l].gene_count,1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);
        fprintf(fp,"%s,%lf,%d,%lf\n",microarray_name[samples[l].gene_count],1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);
        //			result[l].zscore[freq_counter][phase_counter] = (ccorig-csurr_total)/csurr_total_sum;
        mean_counter++;
        phase_counter++;
      }

      {
        double csurr_total=0;
        double csurr_total_sum=0;
        double ccorig = 0;
        double sample_autocorrelation=0;
        double cosine_autocorrelation=0;
        for(x=0;x<rsn;x++){
          csurr_autocorrelation[x] = csurr[x] = 0;
        }
        for(k=0;k<samples[l].timepoint_number;k++){
          ccorig += (samples[l].sample[k] - samples[l].sample_mean) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          sample_autocorrelation += (samples[l].sample[k] - samples[l].sample_mean) * (samples[l].sample[k] - samples[l].sample_mean);
          cosine_autocorrelation += (test_cosine[counter] - test_cosine_mean[mean_counter]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]test_cosine=%f mean=%f\n",i,j,k,test_cosine[counter] ,test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]     sample=%f sample_mean=%f\n",i,j,k,samples[l].sample[k] , samples[l].sample_mean);
          //printf("value[%d][%d]=%lf\n",l,k,samples[l].sample[k]);
          for(m=0;m<rsn;m++){
            csurr[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
            //fprintf(stderr,"csurr[%d]=%lf\n",m,(samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]));
            csurr_autocorrelation[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]);
            //printf("[%d] sample=%lf sample_mean=%lf test_cosine=%lf cosine_mean=%lf csurr=%lf\n",m,samples[l].random_sample[m][k] ,samples[l].sample_mean,test_cosine[counter],test_cosine_mean[counter],csurr[m]);
          }
          //test_cosine[counter] = cos((M_PI*i*time_points[k])+(((M_PI*i)/100) * (double)j));
          counter++;
        }
        ccorig /= ((double)samples[l].timepoint_number - 1);
        ccorig /= (sqrt((sample_autocorrelation) / ((double)samples[l].timepoint_number - 1)) * sqrt((cosine_autocorrelation) / ((double)samples[l].timepoint_number - 1)));
        for(m=0;m<rsn;m++){
          csurr[m] /= ((double)samples[l].timepoint_number - 1);
          csurr[m] /= (sqrt(csurr_autocorrelation[m] / ((double)samples[l].timepoint_number - 1)) * sqrt(cosine_autocorrelation / ((double)samples[l].timepoint_number - 1)));
          //printf("[%lf][%d]csurr_auto=%lf cosine_auto=%lf ccorig=%lf csurr[%d]=%lf\n",i,j,sqrt(csurr_autocorrelation[m]),sqrt(cosine_autocorrelation),m,ccorig,csurr[m]);
          csurr_total += csurr[m];
        }
        csurr_total /= (double)samples[l].random_sample_number;
        for(m=0;m<rsn;m++){
          csurr_total_sum += (csurr[m] - csurr_total) * (csurr[m] - csurr_total);
        }
        csurr_total_sum /= (double)(samples[l].random_sample_number-1);
        csurr_total_sum = sqrt(csurr_total_sum);
        //			printf("ccorig=%lf csurr_total[%d]=%lf\n",ccorig,j,csurr_total);
        //			fprintf(stderr,"[%d][%lf][%d] test_zscore = %lf\n",samples[l].gene_count,1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);
        fprintf(fp,"%s,%lf,%d,%lf\n",microarray_name[samples[l].gene_count],1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);
        //			result[l].zscore[freq_counter][phase_counter] = (ccorig-csurr_total)/csurr_total_sum;
        mean_counter++;
        phase_counter++;
      }

      for(j=1;j <= 50;j++){
        double csurr_total=0;
        double csurr_total_sum=0;
        double ccorig = 0;
        double sample_autocorrelation=0;
        double cosine_autocorrelation=0;
        for(x=0;x<rsn;x++){
          csurr_autocorrelation[x] = csurr[x] = 0;
        }
        for(k=0;k<samples[l].timepoint_number;k++){
          ccorig += (samples[l].sample[k] - samples[l].sample_mean) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          sample_autocorrelation += (samples[l].sample[k] - samples[l].sample_mean) * (samples[l].sample[k] - samples[l].sample_mean);
          cosine_autocorrelation += (test_cosine[counter] - test_cosine_mean[mean_counter]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]test_cosine=%f cosine_mean=%f sample=%f sample_mean=%f\n",i,j,k,test_cosine[counter] ,test_cosine_mean[mean_counter],samples[l].sample[k],samples[l].sample_mean);
          //fprintf(stderr,"[%lf][%d][%d]test_cosine=%f cosine_mean=%f\n",i,j,k,test_cosine[counter] ,test_cosine_mean[mean_counter]);
          //fprintf(stderr,"[%lf][%d][%d]     sample=%f sample_mean=%f\n",i,j,k,samples[l].sample[k] , samples[l].sample_mean);
          //printf("value[%d][%d]=%lf\n",l,k,samples[l].sample[k]);
          for(m=0;m<rsn;m++){
            csurr[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]);
            //fprintf(stderr,"csurr[%d]=%lf\n",m,(samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (test_cosine[counter] - test_cosine_mean[mean_counter]));
            csurr_autocorrelation[m] += (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]) * (samples[l].random_sample[m][k] - samples[l].random_sample_mean[m]);
            //printf("[%d] sample=%lf sample_mean=%lf test_cosine=%lf cosine_mean=%lf csurr=%lf\n",m,samples[l].random_sample[m][k] ,samples[l].sample_mean,test_cosine[counter],test_cosine_mean[counter],csurr[m]);
          }
          //test_cosine[counter] = cos((M_PI*i*time_points[k])+(((M_PI*i)/100) * (double)j));
          counter++;
        }
        ccorig /= ((double)samples[l].timepoint_number - 1);
        ccorig /= (sqrt((sample_autocorrelation) / ((double)samples[l].timepoint_number - 1)) * sqrt((cosine_autocorrelation) / ((double)samples[l].timepoint_number - 1)));
        //			fprintf(stderr,"ccorig=%f\n",ccorig);
        for(m=0;m<rsn;m++){
          csurr[m] /= ((double)samples[l].timepoint_number - 1);
          csurr[m] /= (sqrt(csurr_autocorrelation[m] / ((double)samples[l].timepoint_number - 1)) * sqrt(cosine_autocorrelation / ((double)samples[l].timepoint_number - 1)));
          //fprintf(stderr,"[%lf][%d]csurr_auto=%lf cosine_auto=%lf ccorig=%lf csurr[%d]=%lf\n",i,j,sqrt(csurr_autocorrelation[m]),sqrt(cosine_autocorrelation),m,ccorig,csurr[m]);
          csurr_total += csurr[m];
        }
        csurr_total /= (double)samples[l].random_sample_number;
        for(m=0;m<rsn;m++){
          csurr_total_sum += (csurr[m] - csurr_total) * (csurr[m] - csurr_total);
        }
        csurr_total_sum /= (double)(samples[l].random_sample_number-1);
        csurr_total_sum = sqrt(csurr_total_sum);
        //			fprintf(stderr,"ccorig=%lf csurr_total[%d]=%lf SDsurr=%lf\n",ccorig,j,csurr_total,csurr_total_sum);
        //			fprintf(stderr,"[%d][%lf][%d] test_zscore = %lf\n",samples[l].gene_count,1.0/i,j,(ccorig-csurr_total)/csurr_total_sum);

        //			fprintf(fp,"[%d][%lf][%d] test_zscore = %lf\n",samples[l].gene_count,1.0/i,-j,(ccorig-csurr_total)/csurr_total_sum);
        fprintf(fp,"%s,%lf,%d,%lf\n",microarray_name[samples[l].gene_count],1.0/i,j,(ccorig-csurr_total)/csurr_total_sum);
        //			result[l].zscore[freq_counter][phase_counter] = (ccorig-csurr_total)/csurr_total_sum;
        mean_counter++;
        phase_counter++;
      }
      freq_counter++;
    }
    fclose(fp);
  }
  free(csurr);
}

int main(int argc,char *argv[]){
  double prob_a,prob_b,prob_c,prob_d,prob_e;
  int search_area;
  int vector[5][4];
  double *test_cosine;
  double *test_cosine_mean;
  double *microarray;
  char **microarray_name;
  double ***print_prob;
  double time_points[MAX];
  double sum;
  int a,b,c,d,e;
  int w,x,y,z;
  int i,j,k,l;
  int init_x,init_y,init_w,init_z;
  int init_a;

  int start,end;
  int print_area;
  int start_gene_count,end_gene_count;
  int start_point,end_point;

  // For MPI
  double ****all_print_prob;
  int myid, numprocs;
  int dest,src,tag;
  MPI_Status status;
  MPI_Request request;
  int sample_num,gene_num;
  Sample *samples;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(argc == 3){
	}else{
    exit(1);
  }

  check_expression_file(argv[1],&gene_num,&sample_num);
	global_sigma = atof(argv[2]);
	if(global_sigma < 0){
		global_sigma = SIGMA;
	}

  fprintf(stderr,"%d %d\n",gene_num,sample_num);
  microarray = (double*)malloc(gene_num * sample_num * sizeof(double));
  if(microarray==NULL)exit(1);

  microarray_name = (char **)malloc(gene_num * sizeof(char*));
  if(microarray_name==NULL)exit(1);
  for(i=0;i < gene_num; i++){
    microarray_name[i] = (char*)malloc(MAX * sizeof(char));
    if(microarray_name[i]==NULL)exit(1);
  }
  if(myid == 0)load_expression_file(argv[1],microarray,microarray_name,time_points,gene_num,sample_num);

  fprintf(stderr,"sample_num=%d\n",sample_num);

  test_cosine = (double*)malloc(PHASE_NUM * FREQUENCY_NUM * sample_num * sizeof(double));
  test_cosine_mean = (double*)malloc(PHASE_NUM * FREQUENCY_NUM * sizeof(double));
  fprintf(stderr,"%d\n",PHASE_NUM * FREQUENCY_NUM * sample_num);
  if(test_cosine == NULL || test_cosine_mean == NULL)exit(1);
  for(i=0;i<PHASE_NUM * FREQUENCY_NUM;i++){
    test_cosine_mean[i] = 0;
  }
  if(myid == 0){
    generate_test_cosine(test_cosine,test_cosine_mean,time_points,sample_num);
  }
  // マイクロアレイデータを全体のプロセスにコピー
  fprintf(stderr,"%d\n",myid);
  MPI_Bcast(microarray,gene_num*sample_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
  fprintf(stderr,"[%d]name bcast\n",myid);
  //	MPI_Bcast(microarray_name,gene_num,MPI_CHAR,0,MPI_COMM_WORLD);
  for(i=0;i<gene_num;i++){
    MPI_Bcast(microarray_name[i],MAX,MPI_CHAR,0,MPI_COMM_WORLD);
  }
  MPI_Bcast(test_cosine,PHASE_NUM * FREQUENCY_NUM * sample_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(test_cosine_mean,PHASE_NUM * FREQUENCY_NUM,MPI_DOUBLE,0,MPI_COMM_WORLD);

  if((myid+1) == numprocs){
    start_gene_count = myid * (gene_num / numprocs);
    end_gene_count = gene_num-1;
  }else{
    start_gene_count = myid * (gene_num / numprocs);
    end_gene_count = ((myid+1) * (gene_num / numprocs)) - 1;
  }

  fprintf(stderr,"myid = [%d] start_gene_count=%d end_gene_count=%d\n",myid,start_gene_count,end_gene_count);
  samples = (Sample*)malloc((end_gene_count-start_gene_count+1) * sizeof(Sample));
  for(i=0;i<(end_gene_count-start_gene_count+1);i++){
    int j;
    samples[i].gene_count = start_gene_count+i;
    samples[i].sample = (double*)malloc(sample_num * sizeof(double));
    for(j=0;j<sample_num;j++){
      //			printf("mvalue=%f\n", microarray[(i+start_gene_count)*sample_num+j]);
      samples[i].sample[j] = microarray[(i+start_gene_count)*sample_num+j];
    }
    //		memcpy(samples[i].sample,&microarray[(i+start_gene_count)*sample_num],sample_num*sizeof(double));
    samples[i].average = 0;
    samples[i].sigma = global_sigma;
    samples[i].random_sample_number = 1000;
    samples[i].timepoint_number = sample_num;
    samples[i].random_sample = NULL;
  }
  generate_sample_mean(samples,end_gene_count-start_gene_count+1);
  generate_random_samples(samples,end_gene_count-start_gene_count+1);
  //	calculate_zscore(samples,end_gene_count-start_gene_count+1,test_cosine,test_cosine_mean,microarray_name);
  calculate_minimum_SSR(samples,end_gene_count-start_gene_count+1,test_cosine,test_cosine_mean,microarray_name);

  fprintf(stderr,"End Calculating![%d]\n",myid);
  MPI_Finalize();
}
