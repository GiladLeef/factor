#ifndef FAC_HEADERS
#define FAC_HEADERS

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

// Quadratic sieve integers.
typedef uint32_t uint32; // small size, like a factor base prime number (32-bit)
typedef uint64_t uint64; // medium size, like a factor base prime number squared (64-bit)
typedef int64_t int64; // signed type to perform intermediates computations

// 64-bit factorization structure

typedef struct {
	uint64 prime;
	int power;
} fac64_row;

// Factorization manager's structure

typedef struct {
	struct {
		int tty ;
		int verbose;
		uint64 force;
		uint64 help;
		struct {
			uint64_t seed;
			uint64_t custom;
		} rand;
		uint64 demand[3];
		uint64 timeout;
		const char *input_file;
		const char *output_file;
		char output_format;
		uint64 qs_multiplier;
		uint64 qs_base_size;
		uint64 qs_sieve;
		uint64 qs_error_bits;
		uint64 qs_threshold;
		uint64 qs_alloc_mb;
		uint64 qs_large_prime;
		uint64 qs_poly;
		uint64 qs_laziness;
		//
		uint64 qs_sieve_cutoff; // not documented (see the source code)
		uint64 qs_tick_end; // not documented (see the source code)
	} params;
	uint64 duration_ms ;
	FILE *in;
	FILE *out;
	int code;
	struct {
		size_t row_idx;
		size_t total_rows;
		size_t max_factors;
		size_t max_digits;
		size_t max_bits;
	} scale;
	struct {
		cint num;
		cint tmp[10];
		cint_sheet *sheet;
		uint64 duration_ms;
		uint64_t *seed;
		uint32 trial_start;
		int power;
		const char *input_string;
		char *output_string;
		struct {
			cint num;
			int power;
			int prime;
		} *res;
		struct {
			void *base;
			void *now;
		} mem;
	} session;
} state;

// Quadratic sieve structures

struct qs_relation {
	uint32 id; // smooth relations have a non-zero id.
	cint *X;
	struct {
		uint32 *data;
		uint32 length;
	} Y;
	union {
		struct {
			uint32 *data;
			uint32 length;
		} Z;
		struct qs_relation *next;
	} axis;
	// axis :
	// - "Z" field is used by smooth relations.
	// - "next" is used by partial relations sharing the same variation, grouped within a linked list.
};

typedef struct {

	// computation sheet
	cint_sheet *sheet;

	// numbers that are updated
	struct {
		cint N;
		cint FACTOR;
		cint X;
		cint KEY;
		cint VALUE;
		cint CYCLE;
		cint TEMP[5];
	} vars;

	// polynomial vars
	struct {
		cint A;
		cint B;
		cint C;
		cint D;
		uint32 d_bits;
		uint32 offset;
		struct {
			uint32 x_1;
			uint32 half ;
			uint32 x_2 ;
			uint32 x_3 ;
		} span;
		uint32 gray_max;
		uint32 curves;
	} poly;

	// constants
	struct {
		cint kN;
		cint ONE;
		cint TOO_LARGE_PRIME;
		cint MULTIPLIER;
		cint M_HALF;
	} constants;

	// system
	struct {
		uint32 bytes_allocated;
		void *base;
		void *now;
	} mem;

	// parameters and miscellaneous vars
	struct{
		uint64 start;
		uint64 end;
		uint64 now;
		uint32 tick;
		double prev_pct ;
	} time;
	uint64_t * seed;
	uint64 adjustor;
	uint32 multiplier;
	uint32 n_bits;
	uint32 kn_bits;
	struct {
		uint8_t **positions[2];
		uint8_t *sieve;
		uint8_t *flags;
		uint32 length;
		uint32 length_half;
		uint32 cache_size;
	} m;
	uint32 iterative_list[5];
	uint32 error_bits;
	struct {
		uint32 value;
	} threshold;
	uint32 sieve_again_perms;

	// useful data sharing same length
	struct {
		struct {
			uint32 num;
			uint32 size;
			uint32 sqrt_kN_mod_prime;
			uint32 root[2];
		} *data;
		uint32 largest;
		uint32 length;
	} base;

	// useful data sharing same length
	struct {
		uint32 *A_indexes;
		struct {
			cint B_terms;
			uint32 *A_inv_double_value_B_terms;
			uint32 A_over_prime_mod_prime;
			uint32 prime_index;
			uint64 prime_squared;
		} *data;
		struct {
			uint32 defined;
			uint32 subtract_one;
			uint32 double_value;
		} values;
	} s;

	uint32 *buffer[2]; // proportional to "length of factor base" (medium or large)

	// uniqueness trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager uniqueness[3];

	// data collection made by algorithm
	struct {
		struct qs_relation **data;
		struct {
			uint32 prev;
			uint32 now;
			uint32 peak;
			uint32 needs;
			uint32 capacity;
			uint32 by_partial;
		} length;
		uint64 too_large_prime;
	} relations;

	// pointers to the divisors of N are kept in a flat array
	struct {
		uint32 total_primes;
		uint32 length;
		cint **data;
	} divisors;

	// Lanczos has its own struct
	struct {
		uint32 n_iterations;
		uint32 safe_length;
		struct {
			struct qs_relation *relation;
			uint32 y_length;
		} *snapshot;
	} lanczos;

	state *state;

} qs_sheet;

// Factorization manager, file i/o utilities, misc utilities.
static void printHelp(const char *);
static uint64 getNum(const char *);
static int cliParamMatch(const char *, const char *, const char *);
static int read_key_and_3_values(const char **, state *);
static int read_key_and_2_values(const char **, state *);
static int readKeyValue(const char **, state *);
static int readFlags(const char **, state *);
static void random(cint_sheet *, uint64_t *, cint *, char *, int, int);
static void generateInputFile(state *);
static char *cintString(state *, const cint *);
static inline void inlineCint(cint *, size_t, void **);
static void dupCint(cint *, const cint *, void **);
static inline void intToCint(cint *, uint64);
static uint64 cintToInt(const cint *);
static struct avl_node *avl_cint_inserter(void *, const void *);
static void *memAligned(void *);
static uint64_t getTime(void);
static void debug_print(const state *, int , const char *, ...);
static void display_progress(const char *, double);
static void manager_add_factor(state *, const cint *, int, int);
static void manager_add_simple_factor(state *, uint64, int, int);
static void factorPollardsRho(state *);
static int trialDivide(state *, int, int);
static int anyRootCheck(state *, const cint *, cint *, cint *);
static int perfectPowerCheck(state *, int);
static int primeCheck(state *, int);
static int giveUp(state *, int);
static int factor(state *);
static int validate_input_file(state *);
static size_t prepare_file_descriptors(state *);
static int validate_string_number(const char *, state *);
static void output_csv(state *, int, int, int);
static void output_json(state *, int, int, int);
static void output_default(state *, int, int, int);
static void output(state *, int);
static void prepare_sessions(state *);
static void erase_session(state *);
static void clear_sessions(state *);
static void process_single(state *);
static void process_multi(int, const char **, state *);

// Math functions and system utilities.
static inline uint64 xor_random(uint64 *);
static inline uint64 xor_rand(uint64 *, uint64, uint64);
static int is_prime_4669913(uint32);
static double log_computation(double);
static inline uint32 multiplication_modulo(uint64, uint64, uint32);
static inline uint32 power_modulo(uint64, uint64, uint32);
static int kronecker_symbol(uint32, uint32);
static uint32 tonelli_shanks(uint32, uint32);
static uint32 modular_inverse(uint32, uint32);

// 64-bit factorization.
static int bitSize(uint64);
static uint64 mulMod(uint64, uint64, uint64);
static uint64 powMod(uint64, uint64, uint64);
static int isPrime64bits(uint64);
static uint64 pollardsRho(uint64, uint64_t *);
static uint64 nthRoot(uint64, uint64);
static uint64 squareExtraction(uint64 *, int *);
static void rhoWorker(state *, uint64, fac64_row *);

// Quadratic sieve.
static int quadraticSieve(state *, int);
static int inner_continuation_condition(qs_sheet *);
static int outer_continuation_condition(qs_sheet *);
static void qs_parametrize(qs_sheet *);
static uint32 linear_param_resolution(const double [][2], uint32);
static void qs_initialize_state(qs_sheet *, state *, int);
static void qs_adjust_input_size(qs_sheet *);
static void qs_select_multiplier(qs_sheet *);
static void qs_score_default_multipliers(qs_sheet *, uint32 *, size_t);
static void qs_score_alternative_multipliers(qs_sheet *, uint32 *, size_t);
static void qs_allocate_memory(qs_sheet *);
static void qs_generate_factor_base(qs_sheet *);
static void qs_setup_polynomial_parameters(qs_sheet *);
static void get_started_iteration(qs_sheet *);
static void iteration_part_1(qs_sheet *, const cint *, cint *);
static void iteration_part_2(qs_sheet *, const cint *, cint *);
static void iteration_part_3(qs_sheet *, const cint *, const cint *);
static uint32 iteration_part_4(const qs_sheet *, uint32, uint32 **, cint *);
static void iteration_part_5(qs_sheet *, const cint *, const cint *);
static void iteration_part_6(qs_sheet *, const cint *, const cint *, const cint *, cint *);
static void iteration_part_7(qs_sheet *, uint32, const uint32 *restrict);
static void iteration_part_8(qs_sheet *, uint32, const uint32 *);
static cint * qs_divisors_uniqueness_helper(qs_sheet *, const cint *);
static int qs_register_divisor(qs_sheet *);
static void register_relations(qs_sheet *, const cint *, const cint *, const cint *);
static void register_regular_relation(qs_sheet *, const cint *, const uint32 *const restrict[4]);
static void register_partial_relation(qs_sheet *, const cint *, const cint *, const uint32 *const restrict[4]);
static void qs_factorize_using_null_vectors(qs_sheet *, const uint64_t *restrict);
static int qs_process_remaining_factors(qs_sheet *);

// Linear algebra with Block Lanczos algorithm.
static void lanczos_mul_MxN_Nx64(const qs_sheet *, const uint64_t *, uint32, uint64_t *);
static void lanczos_mul_trans_MxN_Nx64(const qs_sheet *, const uint64_t *, uint64_t *);
static void lanczos_mul_64xN_Nx64(const qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static uint64_t lanczos_find_non_singular_sub(const uint64_t *, const uint64_t *, uint64_t *, uint64_t, uint64_t *);
static void lanczos_mul_Nx64_64x64_acc(qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static void lanczos_mul_64x64_64x64(const uint64_t *, const uint64_t *, uint64_t *);
static void lanczos_transpose_vector(qs_sheet *, const uint64_t *, uint64_t **);
static void lanczos_combine_cols(qs_sheet *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
static void lanczos_build_array(qs_sheet *, uint64_t **, size_t, size_t);
static uint64_t *lanczos_block_worker(qs_sheet *);
static void lanczos_reduce_matrix(qs_sheet *);
static uint64_t *block_lanczos(qs_sheet *);

// Verbose level 0: just factorization, no other messages (default when there is no terminal).
// Verbose level 1: also show task progress in percentage (default when there is a terminal).
// Verbose level 2: also displays final status and duration.
// Verbose level 3: also show maintenance messages.
// Verbose level 4: also show Quadratic Sieve detailed information.

#define DEBUG_LEVEL_2(fmt, ...) debug_print(qs->state, 2, fmt, __VA_ARGS__)
#define DEBUG_LEVEL_3(fmt, ...) debug_print(qs->state, 3, fmt, __VA_ARGS__)
#define DEBUG_LEVEL_4(fmt, ...) debug_print(qs->state, 4, fmt, __VA_ARGS__)

#endif //FAC_HEADERS
