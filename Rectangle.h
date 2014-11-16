#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "Position.h"

struct Rectangle {
	Rectangle() : upperLeft(std::numeric_limits<int8_t>::max(), std::numeric_limits<int8_t>::max()), lowerRight(std::numeric_limits<int8_t>::min(), std::numeric_limits<int8_t>::min()) {}
	Position upperLeft;
	Position lowerRight;
};

#endif
