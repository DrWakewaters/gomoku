#ifndef POSITION_H
#define POSITION_H

struct Position {
	Position() : y(-1), x(-1) {}
	Position(int8_t y, int8_t x) : y(y), x(x) {}
	int8_t y;
	int8_t x;
};

#endif
