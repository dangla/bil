#ifndef NEW_H
#define NEW_H

#include <stdlib.h>

#define New(N,T) \
        ((T*) malloc((N) * sizeof(T)))
        
#endif
