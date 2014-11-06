#ifndef PLY_H
#define PLY_H

#include <vector>
#include "Position.h"

const int16_t NO_VALID_SCORE = std::numeric_limits<int16_t>::max()-1;
const int16_t MAX_SCORE = std::numeric_limits<int16_t>::max();
const int16_t MIN_SCORE = std::numeric_limits<int16_t>::min();

struct Ply {
	Ply() : position(), score(NO_VALID_SCORE) {}
	Ply(Position position, int score) : position(position), score(score) {}
	Position position;
	int16_t score;
};

#endif
