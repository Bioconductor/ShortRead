#ifndef _SHORTREAD_CALL_H_
#define _SHORTREAD_CALL_H_

#ifdef __cplusplus
extern "C" {
#endif

#define NEW_CALL(S, T, NAME, ENV, N)				\
    PROTECT((S) = (T) = allocList((N)));			\
    SET_TYPEOF((T), LANGSXP);                                   \
    SETCAR((T), findFun(install((NAME)), (ENV)));               \
    (T) = CDR((T))
#define CSET_CDR(T, NAME, VALUE)				\
    SETCAR((T), (VALUE));                                       \
    SET_TAG((T), install((NAME)));				\
    (T) = CDR((T))
#define CEVAL_TO(S, ENV, GETS) \
    (GETS) = eval((S), (ENV)); \
    UNPROTECT(1)

#ifdef __cplusplus
}
#endif
#endif                          /* _SHORTREAD_CALL_H_ */
