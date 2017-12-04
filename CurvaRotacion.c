#include <stdio.h>
#include <math.h>
#include <stdlib.h>

float likelihood(float *v, float *v_o, int n){
	float chi;
	float l;
	for(int i = 0; i < n; i++){
		float temp = pow((v_o[i] - v[i])/10.0, 2);
		chi += temp;
	}
	l = exp(-0.5*chi*0.001);
	return l;
}

float * model(float *r_o, float Mb, float Md, float Mh, int n){
	float bb = 0.2497;
	float bd = 5.16;
	float ad = 0.3105;
	float ah = 64.3;
	float *v_m = malloc(n*sizeof(float));
	for(int i = 0; i < n; i++){
		v_m[i] = pow(Mb, 0.5)*r_o[i]/pow(pow(r_o[i], 2) + pow(bb, 2), 0.75) + pow(Md, 0.5)*r_o[i]/pow(pow(r_o[i], 2) + pow(bd + ad, 2), 0.75) + pow(Mh, 0.5)/pow(pow(r_o[i], 2) + pow(ah, 2), 0.25);
	}
	return v_m;
}

float randnormal(float mu){
	
	float x1 = 2*drand48() - 1;
	float x2 = 2*drand48() - 1;
	float x = pow(x1, 2) + pow(x2, 2);
	while(x >= 1 || x == 0){
		x1 = 2*drand48() - 1;
		x2 = 2*drand48() - 1;
		x = pow(x1, 2) + pow(x2, 2);
	}
	
	float ns = x1*pow(-2*log(x)/x, 0.5);
	float n = mu + 10*ns;
	
	return n;
}

int main(){
	FILE *in;
	int lines = 300;
	int iter = 50000;
	float *r_obs = malloc(lines*sizeof(float));
	float *vel_obs = malloc(lines*sizeof(float));
	float *vel_init = malloc(lines*sizeof(float));
	float *vel_prime = malloc(lines*sizeof(float));
	float *best_vel = malloc(lines*sizeof(float));
	float *Mb_walk = malloc(iter*sizeof(float));
	float *Md_walk = malloc(iter*sizeof(float));
	float *Mh_walk = malloc(iter*sizeof(float));
	float *l_walk= malloc(iter*sizeof(float));
	
	/*Cargar datos*/
	char filename1[100] = "RadialVelocities.dat";
	in = fopen(filename1, "r");
	if(!in){
		printf("No se pudo leer el archivo.\n");
		exit(1);
	}
	fscanf(in, "%*[^\n]\n", NULL);
	for(int i = 0; i < lines; i++){
		fscanf(in, "%f %f\n", &r_obs[i], &vel_obs[i]);
	}
	fclose(in);
	
	/*Valores iniciales*/
	Mb_walk[0] = 5000.0;
	Md_walk[0] = 5000.0;
	Mh_walk[0] = 5000.0;
	vel_init = model(r_obs, Mb_walk[0], Md_walk[0], Mh_walk[0], lines);
	l_walk[0] = likelihood(vel_init, vel_obs, lines);
	
	/*Iteracion*/
	for(int i = 1; i< iter; i++){
		float Mb_prime = fabs(randnormal(Mb_walk[i - 1]));
		float Md_prime = fabs(randnormal(Md_walk[i - 1]));
		float Mh_prime = fabs(randnormal(Mh_walk[i - 1]));
		
		vel_init = model(r_obs, Mb_walk[i - 1], Md_walk[i - 1], Mh_walk[i - 1], lines);
    		vel_prime = model(r_obs, Mb_prime, Md_prime, Mh_prime, lines);
		
		float l_prime = likelihood(vel_prime, vel_obs, lines);
		
		float alpha = l_prime/l_walk[i - 1];
		if(alpha >= 1.0){
			Mb_walk[i] = Mb_prime;
			Md_walk[i] = Md_prime;
			Mh_walk[i] = Mh_prime;
			l_walk[i] = l_prime;
		}
		else{
			float beta = drand48();
			if(beta <= alpha){
				Mb_walk[i] = Mb_prime;
				Md_walk[i] = Md_prime;
				Mh_walk[i] = Mh_prime;
				l_walk[i] = l_prime;
			}
			else{
				Mb_walk[i] = Mb_walk[i - 1];
				Md_walk[i] = Md_walk[i - 1];
				Mh_walk[i] = Mh_walk[i - 1];
				l_walk[i] = l_walk[i - 1];
			}
		}
	}
	
	/*Encontramos los parámetros*/
	int index_max = 0;
	float l_max = l_walk[0];
	for(int i = 1; i < iter; i++){
		if( l_walk[i] > l_max){
			l_max = l_walk[i];
			index_max = i;
		}
	}
	
	best_vel = model(r_obs, Mb_walk[index_max], Md_walk[index_max], Mh_walk[index_max], lines);
	
	in = fopen("fit.txt", "w");
	for(int i = 0; i < lines; i++){
		fprintf(in, "%f\n", best_vel[i]);
	}
	
	printf("Mb = %f, Md = %f, Mg = %f\n", Mb_walk[index_max], Md_walk[index_max], Mh_walk[index_max]);
	
	return 0;	
}