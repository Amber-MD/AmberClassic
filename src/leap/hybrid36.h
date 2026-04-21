#ifndef HYBRID36_H
#define HYBRID36_H

#define HY36_WIDTH_4_MIN -999
#define HY36_WIDTH_4_MAX 2436111 /* 10000 + 2*26*36*36*36 - 1 */
#define HY36_WIDTH_5_MIN -9999
#define HY36_WIDTH_5_MAX 87440031 /* 100000 + 2*26*36*36*36*36 - 1 */

/* ------------------------------------------------------------------ */
/* Public error codes                                                   */
/* ------------------------------------------------------------------ */
typedef enum {
    HY36_WRAPPED       = -1, // extra code
    HY36_OK            = 0,
    HY36_OUT_OF_RANGE,
    HY36_INVALID_LITERAL,
    HY36_UNSUPPORTED_WIDTH,
    HY36_INTERNAL_ERROR,
} Hy36Error;

Hy36Error hy36encode(unsigned width, int value, char* result);

Hy36Error hy36decode(unsigned width, const char* s, int* result);

#endif /* HYBRID_36_H */
