#ifndef HASH_H
#define HASH_H

#include <limits>
#include <cstdint>
#include <vector>
#include "Ply.h"

// If we are at a position that is represented by a Hash object hash, the ply bestNextPosition will be tried first (beginning with the first element). After that, other possible plies will be tried.
// The size of a hash object is 2*3 + 2 + 2*1 + 1 + 1 = 12 bytes (trimming a few more bytes is possible, but makes things uglier). The full hash value is 64 bits, but if we make sure to have a hash table of at least 1 << 16 elements, using the least significant 16 bits as index, then we only have to store the 16*3 most significant bits in each hash object.
struct Hash {
	Hash() : hashValue{0, 0, 0}, score(0), bestNextPosition(), bestNextPositionLevel(0), hasValidScore(false) {}
	uint16_t hashValue[3];
	int16_t score;
	Position bestNextPosition;
	uint8_t bestNextPositionLevel;
	bool hasValidScore;
};

#endif
