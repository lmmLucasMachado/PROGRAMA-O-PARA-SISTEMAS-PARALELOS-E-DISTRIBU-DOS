// Define CUDA function directives as spaces
#define __global__
#define __device__
#define _POSIX_C_SOURCE 199309L

#include "prho_vow.h"
#include "hash3_code.h"
#include <omp.h>
#include <time.h>

// Defines that the user can control
#define MAX_RUNS                   100
#define MODULO_NUM_BITS             16
#define NUM_THREADS                  8
#define MAX_EMPTY_SLOT_NOT_FOUND   100

// GLOBAL VARIABLES

// CONSTANT VALUED VARIABLES
const int algorithm = VOW;
const int L = NUM_IT_FUNCTION_INTERVALS;
const int nbits = MODULO_NUM_BITS;
const int nworkers = NUM_THREADS;

// NON-CONSTANT VALUED GLOBAL VARIABLES
long long n_unfruitfulcollisions = 0;
long long kvaluefound;

POINT_T   hashtable[TABLESIZE];
POINT_T   EMPTY_POINT = {0, 0, 0, 0};
POINT_T   Psums[MODULO_NUM_BITS], Qsums[MODULO_NUM_BITS];
POINT_T   X[NUM_THREADS];

IT_POINT_SET itPsetBase;

IT_POINT_SET itPsets[NUM_THREADS];


////////////////////////////////////////////////////////////////////////////
//
// Utility Functions
//
////////////////////////////////////////////////////////////////////////////

void printP(POINT_T X)  {
    printf("(%7lld, %7lld), r = %8lld, s = %8lld", X.x, X.y, X.r, X.s);
}

////////////////////////////////////////////////////////////////////////////
//
// Functions to handle elliptic point aritmetic
//
////////////////////////////////////////////////////////////////////////////

POINT_T add_2P(POINT_T p, long long a, long long mod) {
	POINT_T val;
	long long lambda, num, den, invden;

	if (p.y == 0) {
    	// 2P is equal to O(0, 0)
    	val.x = 0;
    	val.y = 0;
    	return(val);
    }

	num = (3*p.x*p.x + a) % mod;
	while (num < 0) {
		num = num + mod;
	}
	den = 2*p.y;
	while (den < 0) {
        den += mod;
    }
    den = den % mod;

	invden = multiplicative_inverse(den, mod);
	lambda = (num * invden) % mod;

	val.x = (lambda*lambda - p.x - p.x);
	while (val.x < 0) {
	    val.x += mod;
	}
	val.x = val.x % mod;

	val.y = (lambda*(p.x - val.x) - p.y);
		while (val.y < 0) {
	    val.y += mod;
	}
	val.y = val.y % mod;

	return(val);
}


POINT_T add_PQ(POINT_T p, POINT_T q, long long mod)
{
	POINT_T val;
	long long lambda, num, den, invden;

	// Test if p or q is the point at infinite
	if (p.x == 0  &&  p.y == 0)  return(q);
	if (q.x == 0  &&  q.y == 0)  return(p);

	if (p.x == q.x) {
    	val.x = 0;
    	val.y = 0;
    	return(val);
    }

	num = (q.y - p.y);
	while (num < 0) {
		num = num + mod;
	}
	num = num % mod;

	den = q.x - p.x;
	while (den < 0) {
        den += mod;
    }
    den = den % mod;

	invden = multiplicative_inverse(den, mod);
	lambda = (num*invden) % mod;

	val.x = (lambda*lambda - p.x - q.x);
	while (val.x < 0) {
	    val.x += mod;
	}
	val.x = val.x % mod;

	val.y = (lambda*(p.x - val.x) - p.y);
	while (val.y < 0) {
	    val.y += mod;
	}
	val.y = val.y % mod;

	return(val);
}


POINT_T subpoints(POINT_T p, POINT_T q, long long a, long long mod)
{
	POINT_T val;

	q.y = -q.y;
	while (q.y < 0) {
		q.y += mod;
	}

	if (p.x==q.x && p.y==q.y)  val = add_2P(p, a, mod);
	else val = add_PQ(p, q, mod);

	return(val);
}


POINT_T addpoints(POINT_T P, POINT_T Q, long long a, long long p, long long order) {
	POINT_T val;

	if (P.x==Q.x && P.y==Q.y)  val = add_2P(P, a, p);
	else val = add_PQ(P, Q, p);

	if (order != 0) {
        val.r = (P.r + Q.r) % order;
        val.s = (P.s + Q.s) % order;
	}

	return(val);
}


POINT_T multpoint(long long n, POINT_T P, long long a, long long mod, long long order) {
	long long i=1;
	POINT_T val = P;

	if (n == 1) return(P);
	if (n < 0)  { val.x = -1;  val.y = -1;  return(val); }
	if (n == 0) { val.x = 0;   val.y = 0;   return(val); }

	do {
		val = addpoints(val, P, a, mod, order);
		i++;
	} while (i < n);

	return(val);
}


long long multiplicative_inverse(long long num, long long modulo)
{
	long long val=0;
    long long a, b, q, r, x, y;
	long long rest, r1, r2, x1=1, x2=0, y1=0, y2=1;

	if (num == 1) return(1);

    if (num >= modulo) {   a = num;   b = modulo;   }
    else {   a = modulo;   b = num;   }

    r2 = a;
    r1 = b;

    long long i = 1, gcd;

    do {
    	rest = r2 % r1;

    	if (rest == 0) {

    		if (DEBUG) {
    		    printf("\n\n_______________________________________________________________________________\n\n\n");
    		    printf("WILL BREAK:   rest = (r2 %% r1) = 0      r2=%lld    r1=%lld\n\n", r2, r1);
    		    printf("_______________________________________________________________________________\n\n");
    	    }

    	    break;
    	}

    	q = r2/r1;
    	r = rest;
    	x = x1 - q*x2;
    	y = y1 - q*y2;

    	if (DEBUG) {
            printf("\n\n\na[%lld]=%lld   b[%lld]=%lld   q[%lld]=%lld", i, r2, i, r1, i, q);
		    printf("    r[%lld]=%lld    x[%lld]=%lld    y[%lld]=%lld\n\n", i, r, i, x, i, y);
		    printf("%lld = %lld*(%lld) + %lld\n\n", r2, q, r1, r);
		    printf("%lld = %lld*(%lld) + %lld*(%lld)\n", r, a, x, b, y);
	    }

		/* At his point  r = a*x + b*y  */
    	if (r != (a*x + b*y))  {
    		printf("\nError: r != (a * x   +   b * y)\n\n");
    		printf("%lld  !=  (%lld * %lld  +  %lld * %lld)      q=%lld\n", r, a, x, b, y, q);
    		break;
        }

        r2 = r1;
		r1 = r;

		x1 = x2;
		x2 = x;

		y1 = y2;
		y2 = y;

		i++;
    } while (1);

	if (r == 0) {
	    if (num <= modulo) gcd = num;
		else gcd = modulo;
	}
	else gcd = r;

	if (DEBUG) {
	    printf("\ngcd (%lld, %lld) = %lld", num, modulo, gcd);
	    printf("      x = %lld      y = %lld\n\n", x, y);
	    printf("%lld = %lld*(%lld) + %lld*(%lld)\n\n", gcd, a, x, b, y);
    }

	if (gcd == 1)  {
		if (num >= modulo)     /*  a == num,  check x */
		    if (x > 0) {

		    	printf("The multiplicative inverse of %lld mod %lld is %lld\n\n", a, b, x);
		    	if (DEBUG) {
		            printf("The multiplicative inverse of %lld mod %lld is %lld\n\n", a, b, x);
		            printf("(%lld * %lld) mod %lld = %lld\n\n", a, x, b, a*x % b);
		        }
		        val = x;
		    }
		    else {
		    	if (DEBUG)  {
		            printf("The multiplicative inverse of %lld mod %lld is (%lld - %lld) = %lld\n\n", a, b, b, -x, b+x);
		            printf("(%lld * %lld) mod %lld = %lld\n\n", a, b+x, b, (a*(b+x)) % b);
		        }
		        val = b+x;
	        }
	    else     /* b == num,  check y */
	        if (y > 0) {
	        	if (DEBUG)  {
		            printf("The multiplicative inverse of %lld mod %lld is %lld\n\n", b, a, y);
		            printf("(%lld * %lld) mod %lld = %lld\n\n", b, y, a, b*y % a);
		        }
		        val = y;
		    }
		    else {
		    	if (DEBUG) {
		            printf("The multiplicative inverse of %lld mod %lld is (%lld - %lld) = %lld\n\n", b, a, a, -y, a+y);
		            printf("(%lld * %lld) mod %lld = %lld\n\n", b, a+y, a, (b*(a+y)) % a);
		        }
		        val = a+y;
	        }
	}
	else {
		printf("There is no multiplicative inverse of %lld mod %lld\n\n", num, modulo);
		printf("gcd (%lld, %lld) = %lld\n\n", a, b, gcd);
		printf("gcd != 1 --> %lld and %lld are not relatively prime\n\n", a, b);
	}

	return(val);
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// Functions for Pollard Rho & Van Orchoost & Wiener algorithms
//
///////////////////////////////////////////////////////////////////////////////////////////


long long get_group(POINT_T R, int L)  {
    return R.x % L;
}


POINT_T nextpoint(POINT_T P, long long a, long long p, long long order, int L, POINT_T *Triplet) {
    POINT_T Y;
    long long g;

    g = get_group(P, L);
    Y = addpoints(P, Triplet[g], a, p, order);

    return(Y);
}


long long int get_k(long long c1, long long d1, long long c2, long long d2, long long order)
{
    long long  num, den, invden, k;

    num = (c1 - c2);
    while (num < 0)  num = num+order;
    num = num % order;

    den = (d2 - d1);
    while (den < 0)  den = den+order;

    invden = multiplicative_inverse(den, order);
    k = (num*invden) % order;

    return(k);
}

long long get_rand(long long order) {
    long long num = rand() % order;
    return(num);
}

void calc_P_sums(POINT_T P, long long a, long long p, long long order, int nbits) {
    long long i;
    POINT_T Ps=P;

    Psums[0] = Ps;
    for (i=1; i < nbits; i++) {
        Ps = addpoints(Ps, Ps, a, p, 0);
        Psums[i] = Ps;
    }
    return;
}

void calc_Q_sums(POINT_T Q, long long a, long long p, long long order, int nbits) {
    long long i;
    POINT_T Qs=Q;

    Qsums[0] = Qs;
    for (i=1; i < nbits; i++) {
        Qs = addpoints(Qs, Qs, a, p, 0);
        Qsums[i] = Qs;
    }

    return;
}

void calc_Q(POINT_T *pQ, POINT_T *Psums, long long k, int nbits, long long a, long long p) {
    unsigned i, bitmask = 1;

    *pQ = EMPTY_POINT;

    for (i=0; i < nbits; i++) {
        if (k & bitmask) {
            *pQ = addpoints(*pQ, Psums[i], a, p, 0);
        }
        bitmask = (bitmask << 1);
    }

    return;
}

void
rand_itpset(IT_POINT_SET *itpset, POINT_T *Psums, POINT_T *Qsums,
            long long a, long long p, long long order, int L, int nbits,
            int algorithm, int tid) {

    srand((time(NULL)) ^ tid);

    #pragma omp for schedule(dynamic)
    for (long long i=0; i < L; i++) {
        long long r = get_rand(order);
        long long s = get_rand(order);

        itpset->itpoint[i] = EMPTY_POINT;
        uint32_t bitmask = 1;

        for (long long j=0; j < nbits; j++) {
            if (r & bitmask) {
                itpset->itpoint[i] = addpoints(itpset->itpoint[i], Psums[j], a, p, 0);
            }
            if (s & bitmask) {
                itpset->itpoint[i] = addpoints(itpset->itpoint[i], Qsums[j], a, p, 0);
            }
            bitmask = (bitmask << 1);
        }
        itpset->itpoint[i].r = r;
        itpset->itpoint[i].s = s;
    }
    
    return;
}


long long equal(POINT_T X1, POINT_T X2)  {

    if (X1.x != X2.x) return(0);
    if (X1.y != X2.y) return(0);
    if (X1.r != X2.r) return(0);
    if (X1.s != X2.s) return(0);

    return(1);
}



/////////////////////////////////////////////////////////////////////////////
//
// This function attempts to store a new point in the hashtable.
//
/////////////////////////////////////////////////////////////////////////////

void check_slot_and_store(POINT_T X, POINT_T *Y, int token, int *retval)
{
    uint32_t index1 = (uint32_t) X.y, index2 = (uint32_t) ((X.y + 2) >> 1);
    POINT_T H;

    // Modify index2 according to the token value
    if (token > 0) index2 = (index2 + token*index1) % TABLESIZE;

    hashword2((const uint32_t *) &X.x, 1, &index1, &index2);
    index1 = index1 % TABLESIZE;
    index2 = index2 % TABLESIZE;

    // Check slot1; if slot1 is empty, store the point and return
    if (equal(hashtable[index1], EMPTY_POINT)) {
        hashtable[index1] = X;
        *retval = STORED;
        return;
    }
    else {
        // Get the point from slot1
        H = hashtable[index1];

        // Check the point from slot1
        if ((H.x == X.x) && (H.y == X.y))  {
            // If point is the same with different coefficients, return point
            if ((H.r != X.r) || (H.s != X.s)) {
                *Y = H;
                *retval = GOOD_COLLISION;
                return;
            }
        }
        // Point in slot1 was either different or the same with the same coefficients.
        // Check slot2; if slot2 is empty, store the point and return
        if (equal(hashtable[index2], EMPTY_POINT)) {
            hashtable[index2] = X;
            *retval = STORED;
            return;
        }
        else {
            // Get the point from slot2
            H = hashtable[index2];

            // Test the point from slot2
            if ((H.x == X.x) && (H.y == X.y))  {
                // If point in slot 2 is the same with different coefficients, return point
                if ((H.r != X.r) || (H.s != X.s)) {
                    *Y = H;
                    *retval = GOOD_COLLISION;
                    return;
                }
                else
                    // If point is the same with same coefficients report unfruitful
                    *retval = UNFRUITFUL_COLLISION;
                    return;
            }
            else {
                // If point is different report empty slot not found
                *Y = H;
                *retval = EMPTY_SLOT_NOT_FOUND;
                return;
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////////////////
//
// This function calculates a new iteration point set NIPS (composed of L points) by
// multiplication a given iteration point set IPS by an integer mult:
//
//                             NIPS = mult * itSetBase
//
///////////////////////////////////////////////////////////////////////////////////////////

void
mult_PQpset2(IT_POINT_SET *nitpset, IT_POINT_SET *ipset, long long mult, POINT_T *Psums,
             POINT_T *Qsums,int L, int nbits, long long a, long long p, long long order)
{
    int i, j;
    long long r, s;
    uint32_t bitmask;

    for (i=0; i < L; i++) {
        bitmask = 1;
        r = (mult * ipset->itpoint[i].r) % order;
        s = (mult * ipset->itpoint[i].s) % order;

        nitpset->itpoint[i] = EMPTY_POINT;
        for (j=0; j < nbits; j++) {
            if (r & bitmask) {
                nitpset->itpoint[i] = addpoints(nitpset->itpoint[i], Psums[j], a, p, 0);
            }
            if (s & bitmask) {
                nitpset->itpoint[i] = addpoints(nitpset->itpoint[i], Qsums[j], a, p, 0);
            }
            bitmask = (bitmask << 1);
        }
        nitpset->itpoint[i].r = r;
        nitpset->itpoint[i].s = s;
    }

    return;
}


///////////////////////////////////////////////////////////////////////////////////////////
//
// This function calculates the iteration point set IPS (composed of L points) for
// a worker by multiplying the iteration point set base (itSetBase) by an integer m:
//
//                             IPS = m * itSetBase
//
///////////////////////////////////////////////////////////////////////////////////////////

void calc_iteration_point_set2(IT_POINT_SET *itpset, IT_POINT_SET *itpsetbase, POINT_T *Psums, POINT_T *Qsums, int tid,
                              long long a, long long p, long long order, int L, int nbits, int algorithm)

{
    int m;

    if (algorithm == VOW) {
        // Each and every worker has its iteration point set equal to itPsetBase
       *itpset = *itpsetbase;
    }

    else if (algorithm == TOHA) {
        if ((tid % 2) == 0) {
            *itpset = *itpsetbase;
        }
        else {
            // Each odd worker's iteration point set is a multiple of itPsetBase (itPsetBase * (tid+3)/2)
            m = (tid+3)/2;
            mult_PQpset2(itpset, itpsetbase, m, Psums, Qsums, L, nbits, a, p, order);
        }
    }

    else if (algorithm == HARES) {
        // Each and every worker's iteration point set is a multiple of itPsetBase (itPsetBase * (tid+1)
        m = tid+1;
        mult_PQpset2(itpset, itpsetbase, m, Psums, Qsums, L, nbits, a, p, order);
    }

    return;
}




///////////////////////////////////////////////////////////////////////////////////////////
//
// This function calculates the iteration point set IPS (composed of L points) for
// a worker by multiplying the iteration point set base (itSetBase) by an integer m:
//
//                             IPS = m * itSetBase
//
///////////////////////////////////////////////////////////////////////////////////////////

void
calc_iteration_point_set3(IT_POINT_SET *itpset, IT_POINT_SET *itpsetbase,
                          POINT_T *Psums, POINT_T *Qsums, int tid, long long a,
                          long long p, long long order, int L, int nbits,
                          int algorithm)
{
    if (algorithm == VOW) {
        // Each and every worker has its iteration point set equal to itPsetBase
       *itpset = *itpsetbase;
    }

    else if (algorithm == TOHA) {
        if ((tid % 2) == 0) {
            *itpset = *itpsetbase;
        }
        else {
            // Each odd worker's iteration point set is randomly calculated
            // rand_itpset(itpset, Psums, Qsums, a, p, order, L, nbits, algorithm);
        }
    }

    else if (algorithm == HARES) {
        // Each and every worker's iteration point set is randomly calculated
        // rand_itpset(itpset, Psums, Qsums, a, p, order, L, nbits, algorithm);
    }

    return;
}



///////////////////////////////////////////////////////////////////////////////////////////
//
// This function calculates the initial point IP for a worker based on points P and Q:
//
//                                 IP = rP + sQ
//
///////////////////////////////////////////////////////////////////////////////////////////

void initial_point(POINT_T *init_point, POINT_T *Psums, POINT_T *Qsums, int nbits,
                   int gtid, int alg, long long nworkers, long long a,
                   long long p, long long order)
{
    long long i, r, s;
    uint32_t bitmask = 1;
    POINT_T ipoint = { 0, 0, 0, 0 };

    //if ((alg == VOW) || ((alg == TOHA) && (gtid % 2) == 0) || (alg == HARES)) {
    if ((alg == VOW) || (alg == TOHA) || (alg == HARES)) {
        r = get_rand(order);
        s = get_rand(order);

        for (i=0; i < nbits; i++) {
            if (r & bitmask) {
                ipoint = addpoints(ipoint, Psums[i], a, p, 0);
            }
            if (s & bitmask) {
                ipoint = addpoints(ipoint, Qsums[i], a, p, 0);
            }
            bitmask = (bitmask << 1);
        }

        ipoint.r = r;
        ipoint.s = s;

        *init_point = ipoint;
    }
    //else {
    //    *init_point = *(init_point + 1);
    //}

    return;
}



////////////////////////////////////////////////////////////////////////////
//
// This function executes the setup for one worker
//
////////////////////////////////////////////////////////////////////////////

void
setup_worker(long long a, long long p, long long order, const int L, const int nbits, const int id, const int alg) {

    // Calculate the iteration point set for this worker
    calc_iteration_point_set3(&itPsets[id], &itPsetBase, Psums, Qsums, id, a, p, order, L, nbits, alg);

    // Calculate the initial point for this worker
    initial_point(&X[id], Psums, Qsums, nbits, id, alg, nworkers, a, p, order);

    return;
}

/////////////////////////////////////////////////////////////////////////////
//
// This function stores a new point in the hashtable and looks for
// a "good" collision: two points (P1, P2) with equal coordinates
// (X,Y) but with different coefficients (r,s), such that:
//
//      P1 = r1*P + s1*Q,       P2 = r2*P + s2Q,       P1 = P2
//
/////////////////////////////////////////////////////////////////////////////

void handle_newpoint(const POINT_T *X, long long order, long long *key)
{
    int retval, num_emptyslotnotfound = 0;
    POINT_T H;

    do {
        check_slot_and_store(*X, &H, num_emptyslotnotfound, &retval);
        switch (retval)  {
            case GOOD_COLLISION:
                *key = get_k(X->r, X->s, H.r, H.s, order);
                break;
            case EMPTY_SLOT_NOT_FOUND:
                num_emptyslotnotfound++;
                if (num_emptyslotnotfound == MAX_EMPTY_SLOT_NOT_FOUND) {
                    printf("ATTENTION!!! TOO MANY FAILED STORAGE ATTEMPTS -- exiting ...\n");
                    printf("               H = "); printP(H); printf("\n\n");
                    exit(-3);
                }
            case UNFRUITFUL_COLLISION:
                // printf("   --  Unfruitful collision -- continuing ...\n");
                printf("@");
                break;
            case STORED:
                break;
        }

        if ((retval == GOOD_COLLISION) || (retval == UNFRUITFUL_COLLISION)
                                                     || (retval == STORED))  break;

    } while (num_emptyslotnotfound < MAX_EMPTY_SLOT_NOT_FOUND);

    if (retval != GOOD_COLLISION) *key = NO_KEY_FOUND;

    return;
}

////////////////////////////////////////////////////////////////////////////
//
// This function executes one iteration for one worker
//
////////////////////////////////////////////////////////////////////////////

void
worker_it_task(long long a, long long p, long long order, const int L, int id, long long *key) {
    // Calculate the next point
    X[id] = nextpoint(X[id], a, p, order, L, &itPsets[id]);

    // Handle the new point checking if the key has been found
    handle_newpoint(&X[id], order, key);
}

int main() {
    omp_set_num_threads(NUM_THREADS);
    struct timespec start_program, end_program;
    clock_gettime(CLOCK_MONOTONIC, &start_program);
	long long a = 1, b = 44, p, maxorder, k, key;
	int it_number, minits, maxits = 0, i;
	POINT_T Q, P;
    struct timespec start, end;
	double convergence_time, setup_time, total_it_num = 0.0, total_it_time = 0.0, total_su_time = 0.0;
    char *stepdef;

    
	switch (algorithm) {
        case VOW:   printf("\n\nPOLLARD RHO ALGORITHM  --  VOW\n"); stepdef = "all-equal"; break;
        case TOHA:  printf("\n\nPOLLARD RHO ALGORITHM  --  Tortoises and Hares\n"); stepdef = "1/2 equal + 1/2 varying";  break;
        case HARES: printf("\n\nPOLLARD RHO ALGORITHM  --  Tortoises\n"); stepdef = "all-varying"; break;
	}
    

	/////////////////////////////////////////////////////////////////////////
	//
	//     ELLIPTIC CURVE EQUATION (WEIERSTRASS FORM)
	//
	//                 Y^2 = X^3 + aX + b
	//
	/////////////////////////////////////////////////////////////////////////


    printf("\nNumber of bits of the EC prime field (16/20/24/28/32): %d\n\n", nbits);

    //Pick elliptic curve parameters for chosen number of bits of the prime field module p
    switch (nbits)  {
        case 16:
            p = 16747;                         // 16 bits
            maxorder = 16931;        P.x = 1; P.y = 5626;
            k = 8047;
            break;
        case 20:
            p = 1048507;                       // 20 bits
            maxorder = 1049101;      P.x = 373173; P.y = 395411;
            k = 3051;
            break;
        case 24:
            p = 16774421;                       // 24 bits
            maxorder = 16770883;     P.x = 4530807; P.y = 1256865;
            k = 2349;
            break;
        case 28:
            p = 268434997;                      // 28 bits
            maxorder = 268446727;    P.x = 47793986; P.y = 101136283;
            k = 7514;
            break;
       case 32:
            p = 4294966981;                      // 32 bits
            maxorder = 4295084473;   P.x = 3437484969; P.y = 579918983;
            k = 2037;
            break;
        default:
            printf("\n\nnbits is not in the defined set (16, 20, 24, 28, 32)\n\n");
            exit(-1);
    }


	if (a == 1) printf("\nElliptic curve:   y^2 = x^3 + x + %lld      (mod %lld, %d bits)\n\n", b, p, nbits);
	else printf("\nElliptic curve:   y^2 = x^3 + %lldx + %lld      (mod %lld, %d bits)\n\n", a, b, p, nbits);

	printf("Order of the curve = %lld\n\n", maxorder);

	printf("Number of workers = %d\n\n", nworkers);

    // Calculating point Q = kP
    printf("Calculating point Q = kP      (k = %lld)\n", k);
    calc_P_sums(P, a, p, maxorder, nbits);
    calc_Q(&Q, Psums, k, nbits, a, p);
    calc_Q_sums(Q, a, p, maxorder, nbits);

    printf("\n");
    printf("P = (x = %8lld, y = %8lld) (r = %8lld, s = %8lld)         (Base Point)\n", P.x, P.y, P.r, P.s);
    printf("Q = (x = %8lld, y = %8lld) (r = %8lld, s = %8lld)         (Q = kP    )\n\n\n", Q.x, Q.y, Q.r, Q.s);

    // Initialize minits
    minits = maxorder;

    // Loop to calculate the desired secret key (k) MAX_RUNS times, each one
    // with randomly chosen step and starting points for each thread (worker).
    // Remember the number of iterations each run took to find the key (k).
    // Then calculate the expected (average) number of iterations that this
    // calculation takes.
    #pragma omp parallel
    {
        int run, thread_num = omp_get_thread_num();

        // Initialize the random generator
        srand((time(NULL)) ^ thread_num);

        for(run = 1; run <= MAX_RUNS; ++run) {
            #pragma omp single
            {
                // Initialize the number of iterations
                it_number = 1;

                // Start counting the setup time
                clock_gettime(CLOCK_MONOTONIC, &start);
            }

            // Set up the running environment for the search for all workers
            //printf("Run[%3d] setup:    ", run);

            // Calculate the iteration point set base (randomly)
            rand_itpset(&itPsetBase, Psums, Qsums, a, p, maxorder, L, nbits, algorithm, thread_num);
            setup_worker(a, p, maxorder, L, nbits, thread_num, algorithm);
            //printf("   %3d", id+1);
            //fflush(stdout);

            #pragma omp single
            {
                //printf("\nRun[%3d] iterations: ", run);

                // Stop counting the setup time and calculate it
                clock_gettime(CLOCK_MONOTONIC, &end);

                //setup_time = (double)((end - start) / CLOCKS_PER_SEC);
                setup_time = (end.tv_sec - start.tv_sec);
                setup_time += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
                setup_time *= 1000;
                key = NO_KEY_FOUND;

                // Start counting the execution time
                clock_gettime(CLOCK_MONOTONIC, &start);
            }

            #pragma omp reduction(+:it_number)
            {
                for ( ; ; ) {
                    long long keyTemp;
                    worker_it_task(a, p, maxorder, L, thread_num, &keyTemp);
                    if (keyTemp != NO_KEY_FOUND) {
                        key = keyTemp;
                        break;
                    }

                    if (key != NO_KEY_FOUND) {
                        //printf("  (%d its)\n", it_number);
                        break;
                    }
                    it_number = it_number+1;
                }
            }

            #pragma omp single nowait
            {
                // Recover the current time and calculate the execution time
                clock_gettime(CLOCK_MONOTONIC, &end);
                convergence_time = (end.tv_sec - start.tv_sec);
                convergence_time += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
                convergence_time *= 1000;

                // Run converged. Print information for it.
                printf("\nRun[%3d] converged after %4d iterations: k = %lld, nworkers = %4d,\n",
                       run, it_number, key, nworkers);
                printf("         setup time = %6.1lf ms, conv time = %6.1lf ms (itfs \"%s\")\n\n",
                       setup_time, convergence_time, stepdef);

                // Keep track of the minimum number of iterations needed to converge in all runs.
                if (it_number < minits) {
                    minits = it_number;
                }

                if (it_number > maxits) {
                    maxits = it_number;
                }

                total_it_num  = total_it_num  + it_number;
                total_su_time = total_su_time + setup_time;
                total_it_time = total_it_time + convergence_time;
            }

            // Cleanup the hash table
            #pragma omp for nowait
            for (i=0; i < TABLESIZE; i++) {
                hashtable[i] = EMPTY_POINT;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end_program);
    double total_program = (end_program.tv_sec - start_program.tv_sec);
    total_program += (end_program.tv_nsec - start_program.tv_nsec) / 1000000000.0;
    total_program *= 1000;

     // Final statistics
    printf("\n\nFINAL STATISTICS\n\n");
    printf("Runs = %d, Min Its # = %d,  Max Its = %d, Av. Iteration # = %.0lf\n", MAX_RUNS, minits, maxits, total_it_num/MAX_RUNS);
    printf("Av. Setup Time = %.1lf ms, Total Setup Time = %.1lf ms\n", total_su_time/MAX_RUNS, total_su_time);
    printf("Av. Iteration Time = %.1lf ms, Total Iteration Time = %.1lf ms\n", total_it_time/MAX_RUNS, total_it_time);
    printf("Total Time (Setup + Convergence) = %.1lf ms\n", total_su_time + total_it_time);

    printf("\nTempo Total de Execucao do Programa = %.1lf ms\n", total_program);

    fflush(stdout);

	return 0;
}
