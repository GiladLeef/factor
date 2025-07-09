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
typedef uint8_t uint8; // tiny size, 8 bits unsigned
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
		uint64 QSmultiplier;
		uint64 QSbase_size;
		uint64 QSsieve;
		uint64 QSerrorBits;
		uint64 QSthreshold;
		uint64 QSalloc_mb;
		uint64 QSlarge_prime;
		uint64 QSpoly;
		uint64 QSlaziness;
		uint64 QSsieve_cutoff; 
		uint64 QStick_end; 
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

struct QSRelation {
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
		struct QSRelation *next;
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
		cint n;
		cint factor;
		cint x;
		cint key;
		cint value;
		cint cycle;
		cint temp[5];
	} vars;

	// polynomial vars
	struct {
		cint A;
		cint B;
		cint C;
		cint D;
		uint32 dBits;
		uint32 offset;
		struct {
			uint32 x1;
			uint32 half ;
			uint32 x2;
			uint32 x3;
		} span;
		uint32 grayMax;
		uint32 curves;
	} poly;

	// constants
	struct {
		cint kN;
		cint one;
		cint tooLargePrime;
		cint multiplier;
		cint mHalf;
	} constants;

	// system
	struct {
		uint32 bytesAllocated;
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
	uint32 nBits;
	uint32 knBits;
	struct {
		uint8 **positions[2];
		uint8 *sieve;
		uint8 *flags;
		uint32 length;
		uint32 lengthHalf;
		uint32 cacheSize;
	} m;
	uint32 iterativeList[5];
	uint32 errorBits;
	struct {
		uint32 value;
	} threshold;
	uint32 sieveAgainPerms;

	// useful data sharing same length
	struct {
		struct {
			uint32 num;
			uint32 size;
			uint32 sqrtKNModPrime;
			uint32 root[2];
		} *data;
		uint32 largest;
		uint32 length;
	} base;

	// useful data sharing same length
	struct {
		uint32 *aIndexes;
		struct {
			cint bTerms;
			uint32 *aInvDoubleValueBTerms;
			uint32 aOverPrimeModPrime;
			uint32 primeIndex;
			uint64 primeSquared;
		} *data;
		struct {
			uint32 defined;
			uint32 subtractOne;
			uint32 doubleValue;
		} values;
	} s;

	uint32 *buffer[2]; // proportional to "length of factor base" (medium or large)

	// uniqueness trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager uniqueness[3];

	// data collection made by algorithm
	struct {
		struct QSRelation **data;
		struct {
			uint32 prev;
			uint32 now;
			uint32 peak;
			uint32 needs;
			uint32 capacity;
			uint32 byPartial;
		} length;
		uint64 tooLargePrime;
	} relations;

	// pointers to the divisors of N are kept in a flat array
	struct {
		uint32 totalPrimes;
		uint32 length;
		cint **data;
	} divisors;

	// Lanczos has its own struct
	struct {
		uint32 nIterations;
		uint32 safeLength;
		struct {
			struct QSRelation *relation;
			uint32 yLength;
		} *snapshot;
	} lanczos;

	state *state;

} QSSheet;

// Factorization manager, file i/o utilities, misc utilities.
static void printHelp(const char *);
static uint64 getNum(const char *);
static int cliParamMatch(const char *, const char *, const char *);
static int readKeyAnd3Values(const char **, state *);
static int readKeyAnd2Values(const char **, state *);
static int readKeyValue(const char **, state *);
static int readFlags(const char **, state *);
static void rand_cint(cint_sheet *, uint64_t *, cint *, char *, int, int);
static void generateInputFile(state *);
static char *cintString(state *, const cint *);
static inline void inlineCint(cint *, size_t, void **);
static void dupCint(cint *, const cint *, void **);
static inline void intToCint(cint *, uint64);
static uint64 cintToInt(const cint *);
static struct avl_node *avl_cint_inserter(void *, const void *);
static void *memAligned(void *);
static uint64_t getTime(void);
static void debugPrint(const state *, int , const char *, ...);
static void displayProgress(const char *, double);
static void managerAddFactor(state *, const cint *, int, int);
static void managerAddSimpleFactor(state *, uint64, int, int);
static void factorPollardsRho(state *);
static int trialDivide(state *, int, int);
static int anyRootCheck(state *, const cint *, cint *, cint *);
static int perfectPowerCheck(state *, int);
static int primeCheck(state *, int);
static int giveUp(state *, int);
static int factor(state *);
static int validateInputFile(state *);
static size_t prepareFileDescriptors(state *);
static int validateStringNumber(const char *, state *);
static void outputCsv(state *, int, int, int);
static void outputJson(state *, int, int, int);
static void outputDefault(state *, int, int, int);
static void output(state *, int);
static void prepareSessions(state *);
static void eraseSession(state *);
static void clearSessions(state *);
static void processSingle(state *);
static void processMulti(int, const char **, state *);

// Math functions and system utilities.
static inline uint64 xorRandom(uint64 *);
static inline uint64 xorRandint(uint64 *, uint64, uint64);
static int isTinyPrime(uint32);
static double log(double);
static inline uint32 multiplicationModulo(uint64, uint64, uint32);
static inline uint32 powerModulo(uint64, uint64, uint32);
static int kroneckerSymbol(uint32, uint32);
static uint32 tonelliShanks(uint32, uint32);
static uint32 modularInverse(uint32, uint32);

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
static int innerContinuationCondition(QSSheet *);
static int outerContinuationCondition(QSSheet *);
static void parametrize(QSSheet *);
static uint32 linearParamResolution(const double [][2], uint32);
static void initializeState(QSSheet *, state *, int);
static void adjustInputSize(QSSheet *);
static void selectMultiplier(QSSheet *);
static void scoreDefaultMultipliers(QSSheet *, uint32 *, size_t);
static void scoreAlternativeMultipliers(QSSheet *, uint32 *, size_t);
static void allocateMemory(QSSheet *);
static void generateFactorBase(QSSheet *);
static void setupPolynomialParameters(QSSheet *);
static void getStartedIteration(QSSheet *);
static void iterationPart1(QSSheet *, const cint *, cint *);
static void iterationPart2(QSSheet *, const cint *, cint *);
static void iterationPart3(QSSheet *, const cint *, const cint *);
static uint32 iterationPart4(const QSSheet *, uint32, uint32 **, cint *);
static void iterationPart5(QSSheet *, const cint *, const cint *);
static void iterationPart6(QSSheet *, const cint *, const cint *, const cint *, cint *);
static void iterationPart7(QSSheet *, uint32, const uint32 *restrict);
static void iterationPart8(QSSheet *, uint32, const uint32 *);
static cint * divisorsUniquenessHelper(QSSheet *, const cint *);
static int registerDivisor(QSSheet *);
static void registerRelations(QSSheet *, const cint *, const cint *, const cint *);
static void registerRegularRelation(QSSheet *, const cint *, const uint32 *const restrict[4]);
static void registerPartialRelation(QSSheet *, const cint *, const cint *, const uint32 *const restrict[4]);
static void factorizeUsingNullVectors(QSSheet *, const uint64_t *restrict);
static int processRemainingFactors(QSSheet *);

// Linear algebra with Block Lanczos algorithm.
static void mulMxNNx64(const QSSheet *, const uint64_t *, uint32, uint64_t *);
static void mulTransMxNNx64(const QSSheet *, const uint64_t *, uint64_t *);
static void mul64xNNx64(const QSSheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static uint64_t findNonSingularSub(const uint64_t *, const uint64_t *, uint64_t *, uint64_t, uint64_t *);
static void mulNx6464x64Acc(QSSheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static void mul64x6464x64(const uint64_t *, const uint64_t *, uint64_t *);
static void transposeVector(QSSheet *, const uint64_t *, uint64_t **);
static void combineCols(QSSheet *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
static void buildArray(QSSheet *, uint64_t **, size_t, size_t);
static uint64_t *blockWorker(QSSheet *);
static void reduceMatrix(QSSheet *);
static uint64_t *blockLanczos(QSSheet *);

// Verbose level 0: just factorization, no other messages (default when there is no terminal).
// Verbose level 1: also show task progress in percentage (default when there is a terminal).
// Verbose level 2: also displays final status and duration.
// Verbose level 3: also show maintenance messages.
// Verbose level 4: also show Quadratic Sieve detailed information.

#define LOGGER(level, fmt, ...) debugPrint(qs->state, level, fmt, __VA_ARGS__)

#endif //FAC_HEADERS
