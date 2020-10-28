#include "Randomclass.h"

/* This is the initialization routine for the random number generator. NOTE: The seed variables can have values between: 0 <= IJ <= 31328 0 <= KL <= 30081 
The random number sequences created by these two seeds are of sufficient length to complete an entire calculation with. For example, if sveral different groups
are working on different parts of the same calculation, each group could be assigned its own IJ seed. This would leave each group with 30000 choices for the second 
seed. That is to say, this random number generator can create 900 million different subsequences -- with each subsequence having a length of approximately 10^30. */ 
void Random::RandomInitialise(int ij,int kl) 
{ 
	double s,t; int ii,i,j,k,l,jj,m;
	/* Handle the seed range errors First random number seed must be between 0 and 31328 Second seed must have a value between 0 and 30081 */
	if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) { 
		ij = 1802;
		kl = 9373; 
	} 
	i = (ij / 177) % 177 + 2; 
	j = (ij % 177) + 2;
	k = (kl / 169) % 178 + 1; 
	l = (kl % 169);
	for (ii=0; ii<97; ii++) {
		s = 0.0; 
		t = 0.5; 
		for (jj=0; jj<24; jj++) {
			m = (((i * j) % 179) * k) % 179; 
			i = j; 
			j = k;
			k = m;
			l = (53 * l + 1) % 169; 
			if (((l * m % 64)) >= 32)
				s += t;
			t *= 0.5; 
		}
		u[ii] = s; 
	} 
	c = 362436.0 / 16777216.0; 
	cd = 7654321.0 / 16777216.0; 
	cm = 16777213.0 / 16777216.0;
	i97 = 97; 
	j97 = 33; 
	test = TRUE; 
}

/* This is the random number generator proposed by George Marsaglia in Florida State University Report: FSU-SCRI-87-50 */ 
double Random::RandomUniform(void) {
	double uni; /* Make sure the initialisation routine has been called */ 
	if (!test)
		RandomInitialise(1802,9373);
	uni = u[i97-1] - u[j97-1];
	if (uni <= 0.0)
		uni++;
	u[i97-1] = uni; 
	i97--; 
	if (i97 == 0)
		i97 = 97;
	j97--;
	if (j97 == 0) 
		j97 = 97; 
	c -= cd;
	if (c < 0.0)
		c += cm;
	uni -= c;
	if (uni < 0.0)
		uni++;
	return(uni);
} 


/* ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM. THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
  The function returns a normally distributed pseudo-random number with a given mean and standard devaiation. Calls are made to a function subprogram which
  must return independent random numbers uniform in the interval (0,1). The algorithm uses the ratio of uniforms method of A.J. Kinderman and J.F. Monahan
  augmented with quadratic bounding curves. */ 
double Random::RandomGaussian(double mean,double stddev) { 
	double q,u,v,x,y; 
	/* Generate P = (u,v) uniform in rect. enclosing acceptance region Make sure that any random numbers <= 0 are rejected, since 
					  gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0. */ 
	do {
		u = RandomUniform();
		v = RandomUniform();
		if (u <= 0.0 || v <= 0.0)
		{ 
			u = 1.0; 
			v = 1.0; 
		}
		v = 1.7156 * (v - 0.5); /* Evaluate the quadratic form */
		x = u - 0.449871; 
		y = fabs(v) + 0.386595; 
		q = x * x + y * (0.19600 * y - 0.25472 * x); 
		
		/* Accept P if inside inner ellipse */ 
		if (q < 0.27597)
			break; /* Reject P if outside outer ellipse, or outside acceptance region */ 
	}
	while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u)); 
	/* Return ratio of P's coordinates as the normal deviate */ 
	return (mean + stddev * v / u);
}


/* Return random integer within a range, lower -> upper INCLUSIVE */ 
int Random::RandomInt(int lower,int upper) {
	return((int)(RandomUniform() * (upper - lower + 1)) + lower); 
}

/* Return random double within a range, lower -> upper */ 
double Random::RandomDouble(double lower,double upper) { 
	return((upper - lower) * RandomUniform() + lower); 
} 