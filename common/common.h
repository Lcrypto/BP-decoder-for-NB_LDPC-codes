#pragma once

#ifndef LOG_Q
#define LOG_Q					3
#endif
#define Q 						(1<<LOG_Q)

#include "GaloisField.h"

//using namespace std;

#if LOG_Q==2
#define POLYNOM Binary<111>::value
#elif LOG_Q==3
#define POLYNOM Binary<1011>::value
#elif LOG_Q==4
#define POLYNOM Binary<10011>::value
#elif LOG_Q==5
#define POLYNOM Binary<100101>::value
#elif LOG_Q==6
#define POLYNOM Binary<1000011>::value
#elif LOG_Q==8
#define POLYNOM Binary<100011101>::value
#endif

typedef GaloisField<LOG_Q, POLYNOM> Field;
typedef GaloisFieldElement<Field> 	FieldElement;