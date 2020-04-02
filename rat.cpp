//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

#pragma hdrstop

//---------------------------------------------------------------------------
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

using namespace std;
using namespace boost::multiprecision;
using boost::rational;

#define max_n 20

/* A partition is a set of numbers that represent a number of things in each of a set of bins with a specified
 * total across all bins. Unique partitions can be created by sorting the numbers in each partition in a decreasing order.
 * For example the partitions of 5 are 11111,2111,221,311,32,41,5.
 *
 * In this puzzle we use them to represent the number of faces on a die with the same number.
 *
 * The puzzle starts in the 11111... state and is to determine how many rolls of the die it takes to get to the n state.
 * See riddler.
 */

int n_partitions [max_n] [max_n];	/* [i] [j] number of partitions of i things into bins with no more than j things per bin */
double fact_table [max_n];
boost::rational<cpp_int> rat_fact_table [max_n];
int n_comb [max_n] [max_n];			/* [i] [j] is combin (i,j) */
int **partition_list [max_n] [max_n];	/* [i] [j] gives a set of partitions of i things into bins <= j per bin */
int *partition_list_len [max_n] [max_n];	/* [i] [j] [k] gives length of kth partition of i things into bins <= j per bin */

int **perm_table;					/* set of permutations of a partition */

#define n_debugs 10
bool debugs [n_debugs];

#define debug_build_partitions			debugs [0]
#define debug_enumerate_permutations	debugs [1]
#define debug_verbose					debugs [2]
#define debug_trans_prob				debugs [3]
#define debug_solve						debugs [4]

void build_a_bunch_of_tables ()
{   int i;
	int j;
	int k;
	int kk;
	int ipart;
	int isubpart;

	fact_table [0] = 1.0;
	rat_fact_table [0] = 1;
	for (i = 1; i < max_n; i++) {
		fact_table [i] = fact_table [i - 1] * i;
		rat_fact_table [i] = rat_fact_table [i - 1] * i;
	}

	n_comb [0] [0] = 1;
	for (i = 1; i < max_n; i++) {
		n_comb [0] [i] = 0;
	}
	for (i = 1; i < max_n; i++) {
		for (j = 0; j <= i; j++) {
			n_comb [i] [j] = (j > 0 ? n_comb [i - 1] [j - 1] : 0) + n_comb [i - 1] [j];
		}
	}

	for (i = 0; i < max_n; i++) {
		n_partitions [0] [i] = 1;
	}
	for (i = 1; i < max_n; i++) {
		for (j = 1; j < max_n; j++) {
			n_partitions [i] [j] = 0;
			for (k = 1; k <= i && k <= j; k++) {
				n_partitions [i] [j] += n_partitions [i - k] [k];
			}
			if (i == j && debug_build_partitions) {
				printf ("%d %d %d\n", i, j, n_partitions [i] [j]);
			}

		}
	}

	/* Build the partitions. each partition set consists of all integers 1..n followed by all possible partitions n-i, i.
	 * Based on this construction order we can rely on the fact that partition [0] is {111...} and partition [nparts - 1] is {n}.
	*/

	for (i = 0; i < max_n; i++) {
		for (j = 1; j < max_n; j++) {
			if (n_partitions [i] [j] > 0) {
				partition_list [i] [j] = (int **) malloc (n_partitions [i] [j] * sizeof (int *));
				partition_list_len [i] [j] = (int *) malloc (n_partitions [i] [j] * sizeof (int));
				for (k = 0; k < n_partitions [i] [j]; k++) {
					partition_list_len [i] [j] [k] = 0;
					partition_list [i] [j] [k] = (int *) malloc (max_n * sizeof (int));
				  	for (kk = 0; kk < max_n; kk++) {
						partition_list [i] [j] [k] [kk] = 0;
					}
				}
				ipart = 0;
				for (k = 1; k <= i && k <= j; k++) {
				  	for (kk = 0; kk < n_partitions [i - k] [k]; kk++) {
					  	partition_list [i] [j] [ipart] [0] = k;
						partition_list_len [i] [j] [ipart] = 1;
					  	for (isubpart = 0; isubpart < max_n - 1; isubpart++) {
						  	partition_list [i] [j] [ipart] [isubpart + 1] = partition_list [i - k] [k] [kk] [isubpart];
							if (partition_list [i] [j] [ipart] [isubpart + 1] > 0) {
								partition_list_len [i] [j] [ipart] = isubpart + 2;
							}
					  	}
					  	ipart++;
				  	}
				}
				if (debug_build_partitions) {
					printf ("%d %d: ", i, j);
					for (k = 0; k < n_partitions [i] [j]; k++) {
						printf (" [%d]:", partition_list_len [i] [j] [k]);
						for (kk = 0; kk < i; kk++) {
							if (kk > 0)
								printf (",");
							printf ("%d", partition_list [i] [j] [k] [kk]);
						}
					}
					printf ("\n");
				}
			}

		}
	}
}

/* Number of combinations of ntop things into nbot ways. generalized form with a vector of nnbot things.
 * nbot should sum to ntop.
 */

double ncomb (int ntop, int nnbot, int *nbot) {
	int i;
	double r;

	r = fact_table [ntop];
	for (i = 0; i < nnbot; i++) {
		r /= fact_table [nbot [i]];
	}
	return (r);
}

boost::rational<cpp_int> rat_ncomb (int ntop, int nnbot, int *nbot) {
	int i;
	boost::rational<cpp_int> r;

	r = rat_fact_table [ntop];
	for (i = 0; i < nnbot; i++) {
		r /= rat_fact_table [nbot [i]];
	}
	return (r);
}

double my_powi (double x, int n) {
	double r;
	int i;

	r = 1;
	for (i = 0; i < n; i++) {
		r *= x;
	}
	return r;
}

boost::rational<cpp_int> rat_my_powi (boost::rational<cpp_int> x, int n) {
	boost::rational<cpp_int> r;
	int i;

	r = 1;
	for (i = 0; i < n; i++) {
		r *= x;
	}
	return r;
}



/* enumerate all permtations recursively. This is working on the iidx'th number in nums, which has a total of n_nums.
 * Each number in nums is > 0 and in descding order nums [i + 1] <= nums [i], so there may be duplicates.
 * Each permutation must be unique so must avoid producing identical permutations when there are duplicates in nums.
 * We will build each permutation in perm_vec, which has n_out entries.
 * pos_vect has size n_nums and keeps track of where the nums are stored in perm_vec.
 * Perm_vec will be cleared on entry to this, when iidx == 0.
 * We will count the number of permutations in n_perms and save them in perm_set if it is not null.
 */

void enumerate_permutations (int iidx, int n_nums, int *nums, int *pos_vec, int n_out, int *perm_vec, int *n_perms, int **perm_set) {
	int i;
	int iout_start;
	int iout;

	/* check if this is top level entry and clear perm_vec and n_perms */

	if (iidx == 0) {
		for (iout = 0; iout < n_out; iout++) {
			perm_vec [iout] = 0;
		}
		*n_perms = 0;
	}

	if (iidx == n_nums) {
		if (perm_set != NULL) {
			for (i = 0; i < n_out; i++) {
				perm_set [*n_perms] [i] = perm_vec [i];
			}
		}
		(*n_perms)++;
	} else {
		/* Iterate for each possible output position for the number in iidx. To avoid duplicates,
		 * check if the number is the same as previous, and start at the next position after the current pos_vec for it.
		 */

		if (iidx == 0 || nums [iidx] != nums [iidx - 1]) {
			iout_start = 0;
		} else {
			iout_start = pos_vec [iidx - 1] + 1;
		}
		for (iout = iout_start; iout < n_out; iout++) {
			if (perm_vec [iout] == 0) {
				pos_vec [iidx] = iout;
				perm_vec [iout] = nums [iidx];
				enumerate_permutations (iidx + 1, n_nums, nums, pos_vec, n_out, perm_vec, n_perms, perm_set);
				perm_vec [iout] = 0;
			}
		}
	}
}

void test_enum (int nin, int *in_vec, int nout, int **out_vec) {
	int pos_vec [max_n];
	int perm_vec [max_n];
	int nperms;
	int iperm;
	int iout;
	int i;

	nperms = 0;
	enumerate_permutations (0, nin, in_vec, pos_vec, nout, perm_vec, &nperms, out_vec);
	printf ("test perm %d:[", nin);
	for (i = 0; i < nin; i++) {
		if (i > 0) {
			printf (",");
		}
		printf ("%d", in_vec [i]);
	}
	printf ("] to %d has %d perms\n", nout, nperms);
	for (iout = 0; iout < nperms; iout++) {
		for (i = 0; i < nout; i++) {
			printf (" %d", out_vec [iout] [i]);
		}
		printf ("\n");
	}
}


void test_enumerate_perms () {
	int i;
	int in_vec [max_n];
	int **out_vec;

	out_vec = (int **) malloc (10000 * sizeof (int *));
	for (i = 0; i < 10000; i++) {
		out_vec [i] = (int *) malloc (max_n * sizeof (int));
	}

	in_vec [0] = 4;
	in_vec [1] = 2;
	test_enum (2, in_vec, 6, out_vec);

	in_vec [0] = 3;
	in_vec [1] = 2;
	in_vec [2] = 1;
	test_enum (3, in_vec, 6, out_vec);

	in_vec [0] = 2;
	in_vec [1] = 2;
	in_vec [2] = 2;
	test_enum (3, in_vec, 6, out_vec);

}


int count_perms (int n_puz, int plen, int *p) {
	int n_perms = 0;
	int posvec [max_n];
	int permvec [max_n];
	int i;

	enumerate_permutations (0, plen, p, posvec, n_puz, permvec, &n_perms, (int **) NULL);

	return n_perms;
}

/* we need a table to store all permutations of a partition to determine all possible ways
 * one partition can go to another partition. So we might as well determine the largest
 * needed table and alloc it once. Store it as a global since we are hacks.
 */

void init_temp_perm_table (int n_puz, int n_state)
{	int ipart;
	int max_perms;
	int n_perms;

	max_perms = 0;
	for (ipart = 0; ipart < n_state; ipart++) {
		n_perms = count_perms (n_puz, partition_list_len [n_puz] [n_puz] [ipart], partition_list [n_puz] [n_puz] [ipart]);
		if (n_perms > max_perms) {
			max_perms = n_perms;
		}
	}
	if (debug_verbose) {
		printf ("max perms %d\n", max_perms);
	}
	perm_table = (int **) malloc (max_perms * sizeof (int *));
	for (ipart = 0; ipart < max_perms; ipart++) {
		perm_table [ipart] = (int *) malloc (n_puz * sizeof (int));
	}
	printf ("%d max perms\n", max_perms);
}

double compute_trans_prob (int n_puz, int sfrom, int sto) {
	int *pfrom;
	int *pto;
	int npfrom;
	int npto;
	int nperms;
	int iperm;
	int ito;
	int ifrom;
	int i;
	int pos_vec [max_n];
	int perm_vec [max_n];
	double sum_prob; 		/* total prob of the transition */
	double perm_prob;		/* prob of one particular permutation */

	/* find the partitions for the from and to states */

	pfrom = partition_list [n_puz] [n_puz] [sfrom];
	npfrom = partition_list_len [n_puz] [n_puz] [sfrom];
	pto = partition_list [n_puz] [n_puz] [sto];
	npto = partition_list_len [n_puz] [n_puz] [sto];

	if (debug_trans_prob) {
		printf ("trans from");
		for (i = 0; i < npfrom; i++) {
			printf (" %d", pfrom [i]);
		}
		printf (" to");
		for (i = 0; i < npto; i++) {
			printf (" %d", pto [i]);
		}
		printf ("\n");
	}

	/* find all  permutations of the to state */

	enumerate_permutations (0, npto, pto, pos_vec, npfrom, perm_vec, &nperms, perm_table);
	sum_prob = 0;
	for (iperm = 0; iperm < nperms; iperm++) {
		perm_prob = 1.0;
		if (debug_trans_prob) {
			printf ("perm");
			for (i = 0; i < npfrom; i++) {
				printf (" %d", perm_table [iperm] [i]);
			}
		}
		/* compute prob based on number of each entry in the from partition */

		for (ifrom = 0; ifrom < npfrom; ifrom++) {
			perm_prob *= my_powi (pfrom [ifrom] / (double) n_puz,perm_table [iperm] [ifrom]);
		}
		perm_prob *= ncomb (n_puz, npto, pto);
		if (debug_trans_prob) {
			printf (" = %g\n", perm_prob);
		}
		sum_prob += perm_prob;
	}
	if (debug_trans_prob) {
		printf ("trans prob = %g\n", sum_prob);
	}
	return sum_prob;

}

boost::rational<cpp_int> rat_compute_trans_prob (int n_puz, int sfrom, int sto) {
	int *pfrom;
	int *pto;
	int npfrom;
	int npto;
	int nperms;
	int iperm;
	int ito;
	int ifrom;
	int i;
	int pos_vec [max_n];
	int perm_vec [max_n];
	boost::rational<cpp_int> sum_prob; 		/* total prob of the transition */
	boost::rational<cpp_int> perm_prob;		/* prob of one particular permutation */
	boost::rational<cpp_int> ptop;
	boost::rational<cpp_int> pbot;

	/* find the partitions for the from and to states */

	pfrom = partition_list [n_puz] [n_puz] [sfrom];
	npfrom = partition_list_len [n_puz] [n_puz] [sfrom];
	pto = partition_list [n_puz] [n_puz] [sto];
	npto = partition_list_len [n_puz] [n_puz] [sto];

	if (debug_trans_prob) {
		printf ("trans from");
		for (i = 0; i < npfrom; i++) {
			printf (" %d", pfrom [i]);
		}
		printf (" to");
		for (i = 0; i < npto; i++) {
			printf (" %d", pto [i]);
		}
		printf ("\n");
	}

	/* find all  permutations of the to state */

	enumerate_permutations (0, npto, pto, pos_vec, npfrom, perm_vec, &nperms, perm_table);
	sum_prob = 0;
	for (iperm = 0; iperm < nperms; iperm++) {
		perm_prob = 1;
		if (debug_trans_prob) {
			printf ("perm");
			for (i = 0; i < npfrom; i++) {
				printf (" %d", perm_table [iperm] [i]);
			}
		}
		/* compute prob based on number of each entry in the from partition */

		for (ifrom = 0; ifrom < npfrom; ifrom++) {
			ptop = pfrom [ifrom];
			pbot = n_puz;
			perm_prob *= rat_my_powi (ptop / pbot, perm_table [iperm] [ifrom]);
		}
		perm_prob *= rat_ncomb (n_puz, npto, pto);
		if (debug_trans_prob) {
			printf (" = %g\n", perm_prob);
		}
		sum_prob += perm_prob;
	}
	if (debug_trans_prob) {
		cout << "trans prob = " << sum_prob << "\n";
	}
	return sum_prob;

}

void solve_matrix ( int n, double **a, double *rhs, double *sol ) {
	int i, j, k;
	double pivotval;
	double t;
#ifdef findpivot
	int pivotrow;
#endif

	for (i = 0; i < n; i++)
		sol [i] = rhs [i];
	for (i = 0; i < n; i++)
	{
#ifdef findpivot
		pivotrow = i;
		pivotval = a [i] [i];
		for (j = i + 1; j < n; j++)
		{	if (abs (a [j] [i]) > pivotval)
			{	pivotval = my_abs (a [j] [i]);
				pivotrow = j;
			}
		}
		for (j = 0; j < n; j++)
		{	t = a [i] [j];
			a [i] [j] = a [pivotrow] [j];
			a [pivotrow] [j] = t;
		}
		t = sol [i];
		sol [i] = sol [pivotrow];
		sol [pivotrow] = t;
#endif
		pivotval = 1.0 / a [i] [i];
		for (j = 0; j < n; j++)
		{	a [i] [j] *= pivotval;
		}
		sol [i] *= pivotval;
		for (j = i + 1; j < n; j++)
		{	pivotval = a [j ] [i];
			a [j] [i] = 0.0;
			for (k = i + 1; k < n; k++)
			{	a [j] [k] -= pivotval * a [i] [k];
			}
			sol [j] -= pivotval * sol [i];
		}
	}
	for (i = n - 1; i >= 0; i--)
	{	for (j = i - 1; j >= 0; j--)
		{	pivotval = a [j] [i];
			a [j] [i] = 0.0;
			sol [j] -= pivotval * sol [i];
		}
	}
}

void rat_solve_matrix ( int n, boost::rational<cpp_int> **a, boost::rational<cpp_int> *rhs, boost::rational<cpp_int> *sol ) {
	int i, j, k;
	boost::rational<cpp_int> pivotval;

	for (i = 0; i < n; i++)
		sol [i] = rhs [i];
	for (i = 0; i < n; i++)
	{
		pivotval = 1 / a [i] [i];
		for (j = 0; j < n; j++)
		{	a [i] [j] *= pivotval;
		}
		sol [i] *= pivotval;
		for (j = i + 1; j < n; j++)
		{	pivotval = a [j ] [i];
			a [j] [i] = 0;
			for (k = i + 1; k < n; k++)
			{	a [j] [k] -= pivotval * a [i] [k];
			}
			sol [j] -= pivotval * sol [i];
		}
	}
	for (i = n - 1; i >= 0; i--)
	{	for (j = i - 1; j >= 0; j--)
		{	pivotval = a [j] [i];
			a [j] [i] = 0;
			sol [j] -= pivotval * sol [i];
		}
	}
}


void solve_puz (int n_puz, int niters) {
	int i;
	int j;
	int n_state;
	double **state_trans_prob;	/* [i] [j] is prob of transitioning from state i to state j */
	int ifrom;
	int ito;
	double rowsum;
	double *distance;
	double *distance_next;
	double delta0;
	int iiter;
	double *rhs;
	double *sol;

	n_state = n_partitions [n_puz] [n_puz];
	printf ("puz %d has %d states\n", n_puz, n_state);
	if (debug_verbose) {
		printf ("puz %d has %d states\n", n_puz, n_state);
	}
	state_trans_prob = (double **) malloc (n_state * sizeof (double *));
	rhs = (double *) malloc (n_state * sizeof (double));
	sol = (double *) malloc (n_state * sizeof (double));
	for (i = 0; i < n_state; i++) {
		state_trans_prob [i] = (double *) malloc (n_state * sizeof (double));
	}

	init_temp_perm_table (n_puz, n_state);

	for (ifrom = 0; ifrom < n_state; ifrom++) {
		for (ito = 0; ito < n_state; ito++) {
			state_trans_prob [ifrom] [ito] = compute_trans_prob (n_puz, ifrom, ito);
		}
	}
	if (debug_trans_prob) {
		printf ("prob table\n");
		for (ifrom = 0; ifrom < n_state; ifrom++) {
			rowsum = 0;		/* row should sum to 1 or we have fucked up */
			for (ito = 0; ito < n_state; ito++) {
				printf (" %12g", state_trans_prob [ifrom] [ito]);
				rowsum += state_trans_prob [ifrom] [ito];
			}
			printf ("	sum %12g\n", rowsum);
		}
	}

	/* now solve for the distances using iteration because I'm too lazy to do direct solve for now */
	if (niters > 0) {
		distance = (double *) malloc (n_state * sizeof (double));
		distance_next = (double *) malloc (n_state * sizeof (double));

		for (ifrom = 0; ifrom < n_state; ifrom++) {
			distance [i] = 1.0;		/* using 0 would be more realistic but this gives it a bit of head start since it is >= 1 */
		}
		distance [n_state - 1] = 0.0;
		for (iiter = 0; iiter < niters; iiter++) {
			if (debug_solve) {
				printf ("iter %6d", iiter);
			}
			for (ifrom = 0; ifrom < n_state - 1; ifrom++) {
				distance_next [ifrom] = 0;
				for (ito = 0; ito < n_state; ito++) {
					distance_next [ifrom] += state_trans_prob [ifrom] [ito] * distance [ito];
				}
				distance_next [ifrom] += 1;
				if (debug_solve) {
					printf (" %6g", distance_next [ifrom]);
				}
			}
			delta0 = distance_next [0] - distance [0];
			if (debug_solve) {
				printf (" delta0 = %12.6g\n", delta0);
			}
			for (ifrom = 0; ifrom < n_state - 1; ifrom++) {
				distance [ifrom] = distance_next [ifrom];
			}

			if (debug_solve) {
				printf ("\n");
			}

		}
		printf ("answer = %20.13g delta0 = %g\n", distance [0], delta0);
	}

	for (i = 0; i < n_state - 1; i++) {
		state_trans_prob [i] [i] -= 1;
		rhs [i] = -1;
	}
	rhs [n_state - 1] = 0;
	solve_matrix (n_state, state_trans_prob, rhs, sol);
	printf ("answer = %20.13g\n", sol [0]);

}

void rat_solve_puz (int n_puz, int niters) {
	int i;
	int j;
	int n_state;
	boost::rational<cpp_int> **state_trans_prob;	/* [i] [j] is prob of transitioning from state i to state j */
	int ifrom;
	int ito;
	boost::rational<cpp_int> rowsum;
	boost::rational<cpp_int> *distance;
	int iiter;
	boost::rational<cpp_int> *rhs;
	boost::rational<cpp_int> *sol;

	n_state = n_partitions [n_puz] [n_puz];
	printf ("puz %d has %d states\n", n_puz, n_state);
	if (debug_verbose) {
		printf ("puz %d has %d states\n", n_puz, n_state);
	}
	state_trans_prob = new boost::rational<cpp_int> * [n_state];
	rhs = new boost::rational<cpp_int> [n_state];
	sol = new boost::rational<cpp_int> [n_state];
	for (i = 0; i < n_state; i++) {
		state_trans_prob [i] = new boost::rational<cpp_int> [n_state];
	}

	init_temp_perm_table (n_puz, n_state);

	for (ifrom = 0; ifrom < n_state; ifrom++) {
		for (ito = 0; ito < n_state; ito++) {
			state_trans_prob [ifrom] [ito] = rat_compute_trans_prob (n_puz, ifrom, ito);
		}
	}
	if (debug_trans_prob) {
		printf ("prob table\n");
		for (ifrom = 0; ifrom < n_state; ifrom++) {
			rowsum = 0;		/* row should sum to 1 or we have fucked up */
			for (ito = 0; ito < n_state; ito++) {
				cout << " " << state_trans_prob [ifrom] [ito];
				rowsum += state_trans_prob [ifrom] [ito];
			}
			cout << "      sum " << rowsum << "\n";
		}
	}


	for (i = 0; i < n_state - 1; i++) {
		state_trans_prob [i] [i] -= 1;
		rhs [i] = -1;
	}
	rhs [n_state - 1] = 0;
	rat_solve_matrix (n_state, state_trans_prob, rhs, sol);
	cout << "rational answer =  " << sol [0] << "\n";

}

void solve_puz_monte (int npuz, int nmonte) {
	int iiter;
	int nrolls;
	int nrolls_test;
	int die_state [max_n];
	int die_state_next [max_n];
	bool done;
	int i;

	nrolls = 0;
	for (iiter = 0; iiter <nmonte; iiter++) {
		for (i = 0; i < npuz; i++) {
			die_state [i] = i;
		}
		done = false;
		nrolls_test = 0;
		while (!done) {
			for (i = 0; i < npuz; i++) {
				die_state_next [i] = die_state [rand () % npuz];
			}
			nrolls_test++;
			done = true;
			for (i = 0; i < npuz; i++) {
				die_state [i] = die_state_next [i];
				if (die_state [i] != die_state [0]) {
					done = false;
				}
			}
		}
//		printf ("%d\n", nrolls_test);
		nrolls += nrolls_test;
	}
	printf ("monte %g\n", (double) nrolls / nmonte);
}

#pragma argsused

int main(int argc, char* argv[])
{   int i;
	int n_puz;
	int niters = 0;
	int nmonte = 0;

	for (i = 0; i < n_debugs; i++) {
		debugs [i] = false;
	}
	if (argc >= 2 && strcmp (argv [1], "-d") == 0) {
		debug_build_partitions = true;
		debug_enumerate_permutations = true;
		debug_verbose = true;
		debug_trans_prob = true;
		debug_solve = true;
		argc--; argv++;
	}

	if (argc < 2 || sscanf (argv [1], "%d", &n_puz) != 1 || n_puz <= 0 || n_puz >= max_n) {
		fprintf (stderr, "bad arg\n");
		exit (1);
	}
	build_a_bunch_of_tables ();

	if (debug_enumerate_permutations) {
		test_enumerate_perms ();
	}

	solve_puz (n_puz, niters);
	rat_solve_puz (n_puz, niters);

	if (nmonte > 0) {
		solve_puz_monte (n_puz, nmonte);
	}

	return 0;
}
//---------------------------------------------------------------------------
